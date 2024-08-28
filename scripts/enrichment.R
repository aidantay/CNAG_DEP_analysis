#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

library(argparse)
library(topGO)
library(ggplot2)
library(dplyr)

# Rscript scripts/enrichment/main.R -i data/experiment_3/CBD_V_vs_Control.tsv -u data/experiment_3/GO_KEGG/uniprot/geneID2GO_all.tsv -o all_CBD_V_vs_Control -p &
# Rscript scripts/enrichment/main.R -i data/experiment_3/CBD_vs_Control.tsv -u data/experiment_3/GO_KEGG/uniprot/geneID2GO_all.tsv -o all_CBD_vs_Control -p &
# Rscript scripts/enrichment/main.R -i data/experiment_3/CBD_vs_CBD_V.tsv -u data/experiment_3/GO_KEGG/uniprot/geneID2GO_all.tsv -o all_CBD_vs_CBD_V -p &

# Rscript scripts/enrichment/main.R -i data/experiment_3/CBD_V_vs_Control.tsv -u data/experiment_3/GO_KEGG/uniprot/geneID2GO_identified.tsv -o identified_CBD_V_vs_Control -p &
# Rscript scripts/enrichment/main.R -i data/experiment_3/CBD_vs_Control.tsv -u data/experiment_3/GO_KEGG/uniprot/geneID2GO_identified.tsv -o identified_CBD_vs_Control -p &
# Rscript scripts/enrichment/main.R -i data/experiment_3/CBD_vs_CBD_V.tsv -u data/experiment_3/GO_KEGG/uniprot/geneID2GO_identified.tsv -o identified_CBD_vs_CBD_V -p &

# Rscript scripts/enrichment/main.R -i data/experiment_3/CBD_V_vs_Control.tsv -u data/experiment_3/GO_KEGG/uniprot/geneID2GO_quantified.tsv -o quantified_CBD_V_vs_Control -p &
# Rscript scripts/enrichment/main.R -i data/experiment_3/CBD_vs_Control.tsv -u data/experiment_3/GO_KEGG/uniprot/geneID2GO_quantified.tsv -o quantified_CBD_vs_Control -p &
# Rscript scripts/enrichment/main.R -i data/experiment_3/CBD_vs_CBD_V.tsv -u data/experiment_3/GO_KEGG/uniprot/geneID2GO_quantified.tsv -o quantified_CBD_vs_CBD_V -p &

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

.getParser <- function() {
    p <- ArgumentParser()
    p$add_argument("-i", help="ID list or DE results file (.tsv)", type="character", required=TRUE)
    p$add_argument("-u", help="ID->GO mapping file (.tsv)", type="character", required=TRUE)
    p$add_argument("-o", help="Output directory", type="character", required=TRUE)

    ## DE parameters
    p$add_argument("-lfc", help="Log2 fold change threshold (default: 1)", default=1)
    p$add_argument("-pval", help="DE p-value threshold (default: 0.05)", default=0.05)
    p$add_argument("-p", help="Filter on padj (default: True)", action='store_false')
    return (p)
}

.isValidArgs <- function(p, args) {
    if (is.null(args$i) | is.null(args$u) | is.null(args$o)) {
        p$print_help()
        stop("Must supply input file, gene universe and and output directory.", call.=FALSE)
    }
}

.filterDE <- function(df, lfc=1, pval=0.05, withPadj=TRUE) {
    ## With multiple correction, we can be more confident of an effect.
    ## Without multiple correction, we can't say there's no effect. Instead, only a possible effect.
    cond1 <- (abs(df$log2FoldChange) > lfc)
    if (withPadj) {
      cond2 <- (df$padj < pval)

    } else {
      cond2 <- (df$pvalue < pval)
    }

    df <- df[which(cond1 & cond2), ]
    return (df)
}

runGOAnalysis <- function(idsOfInterest, ID2GO) {
    ## Specify which ids in the "id universe" are ids of interest
    ID2GO <- readMappings(file=ID2GO, sep='\t')
    idUniverse <- names(ID2GO)
    idList <- factor(as.integer(idUniverse %in% idsOfInterest))
    names(idList) <- idUniverse

    ## Run GO enrichment
    go_cats <- c('BP', 'CC', 'MF')
    res <- lapply(go_cats, \(x) getTopGOTable(x, idList, ID2GO))
    res <- res[!is.na(res)]
    res <- lapply(res, \(x) insertRatioCols(x, length(idsOfInterest), length(ID2GO)))
    return (res)
}

getTopGOTable <- function(GO_C, idList, ID2GO) {
    print(GO_C)

    ## Construct the TopGO structure
    GOdata <- new(
        "topGOdata", ontology=GO_C,
        allGenes=idList, nodeSize=1,
        annot=annFUN.gene2GO, gene2GO=ID2GO)
    res <- getGOTable(GOdata)
    if (dim(res)[[1]] == 0) {
        return (NA)
    }

    ## Find ids of interest associated with each GO term
    idsOfInterest <- names(idList[idList == 1])
    res <- insertIDCols(GOdata, res, idsOfInterest)

    ## Format table
    res$Category <- GO_C
    res <- (
        res
        %>% relocate(GO.ID)
        %>% relocate(Category)
        %>% relocate(IDs, .after=Term)
        %>% relocate(Significant, .after=IDs)
        %>% rename(Query='Significant', Background='Annotated')
    )

    return (res)
}

getGOTable <- function(GOdata) {
    ## Construct test statistic & find significant GO Terms
    res.classicFisher  <- runTest(GOdata, algorithm='classic', statistic='fisher')
    res.classicKS      <- runTest(GOdata, algorithm='classic', statistic='ks')

    res <- GenTable(
        GOdata,
        classicFisher=res.classicFisher, classicKS=res.classicKS,
        orderBy="classicFisher", ranksOf="classicFisher",
        topNodes=length(usedGO(GOdata)),
        numChar=1000)

    res$classicFisher_padj  <- p.adjust(res$classicFisher, method='fdr')
    res$classicKS_padj      <- p.adjust(res$classicKS, method='fdr')

    res$classicFisher  <- as.numeric(res$classicFisher)
    res$classicKS      <- as.numeric(res$classicKS)

    ## Remove rows that don't contain our ids of interest
    res <- res[res$Significant > 0, ]
    return (res)
}

insertIDCols <- function(GOdata, res, idsOfInterest) {
    ## Find all GO terms associated with the ids of interest
    GO2ID <- genesInTerm(GOdata, res$GO.ID)
    GOOfInterest <- lapply(GO2ID, \(x) intersect(x, idsOfInterest))
    GOOfInterest <- GOOfInterest[lapply(GOOfInterest,length)>0]
    GOOfInterest <- lapply(GOOfInterest, paste, collapse=",")

    ## Construct a table
    GOOfInterest <- as.data.frame(as.matrix(GOOfInterest))
    colnames(GOOfInterest) <- c("IDs")
    GOOfInterest$IDs   <- as.character(GOOfInterest$IDs)
    GOOfInterest$GO.ID <- rownames(GOOfInterest)

    ## Join the GO results with the ids of interest
    res <- merge(res, GOOfInterest, by='GO.ID', all.x=TRUE)
    res <- res[order(res$classicFisher_padj),]
    return (res)
}

insertRatioCols <- function(res, qCount, bgCount) {
    res$QryRatio <- paste(as.character(res$Query), as.character(qCount), sep='/')
    res$BgRatio  <- paste(as.character(res$Background), as.character(bgCount), sep='/')
    res <- (
        res
        %>% relocate(QryRatio, .after=Query)
        %>% relocate(BgRatio, .after=Background)
    )
    res <- subset(res, select=-c(Query, Background))
    return (res)
}

filterGO <- function(df, col='classicFisher_padj') {
    cond1 <- df[[col]] < 0.05
    cond2 <- df$Query >= 0
    df <- df[which(cond1 & cond2), ]
    return (df)
}

getGOPlot <- function(df, col='classicFisher_padj', how='dotplot') {
    filt_de <- df
    filt_de$GeneRatio <- sapply(filt_de$QryRatio, \(x) eval(parse(text=x)))
    filt_de$Query <- lapply(strsplit(filt_de$QryRatio, '/'), \(x) x[1])
    filt_de$Query <- as.numeric(filt_de$Query)
    filt_de <- filterGO(filt_de, col)
    p <- getGODotplot(filt_de, col)
    return(p)
}

getGODotplot <- function(df, col) {
    p <- (
        ggplot(
            df, aes(
                x=Query,
                y=reorder(Term, as.numeric(as.factor(Category)) + .data[[col]])
            )
        )
        + geom_point(aes(color=.data[[col]], size=GeneRatio))
        + geom_vline(xintercept=c(5), linetype='dashed', color='grey') ## Add vertical line
        + scale_color_gradient(low="blue", high="red")
        + scale_size(name="Gene Ratio")
        + geom_tile(aes(width=Inf, fill=Category), alpha=0.15)
        + labs(x="Count", y=NULL, color=col)
    )
    return(p)
}

#------------------- Main -----------------------------------#

# Parse command-line arguments
p    <- .getParser()
args <- p$parse_args()
.isValidArgs(p, args)
print(args)

## Read data
de_df <- read.table(args$i, sep='\t', header=TRUE, quote="")
de_df$UniprotAcc <- as.character(de_df$UniprotAcc)
de_df$text       <- unlist(lapply(strsplit(de_df$id, ' \\| '), \(x) x[2]))
de_df$id <- unlist(lapply(strsplit(de_df$id, '\\|'), \(x) x[1]))

## Specify which genes are +ve/-ve expressed
filt_de <- .filterDE(de_df, args$lfc, args$pval, args$p)

## Run GO Analysis
all_res <- runGOAnalysis(filt_de$UniprotAcc, args$u)
all_res <- bind_rows(all_res)
all_go_plot <- getGOPlot(all_res, col='classicKS_padj')

dir.create(args$o, showWarnings=FALSE, recursive=TRUE)
ggsave(file.path(args$o, 'GO_All.pdf'), all_go_plot, width=20, height=20)
write.table(all_res, file.path(args$o, 'GO_All.tsv'), quote=FALSE, sep='\t', row.names=FALSE)