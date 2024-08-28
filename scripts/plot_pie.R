#------------------- Description & Notes --------------------#

# Rscript plot_pie.R \
#    -i CBD_V_vs_Control.tsv \
#    -g GO_CBD_V_vs_Control.tsv
#    -o out

#------------------- Dependencies ---------------------------#

library(argparse)
library(ggplot2)
library(tidyverse)

#------------------- Constants ------------------------------#

BP <- list(
    'metabolic process'    = 'Metabolic, catabolic and biosynthetic process',
    'catabolic process'    = 'Metabolic, catabolic and biosynthetic process',
    'biosynthetic process' = 'Metabolic, catabolic and biosynthetic process',
    'cellular process'     = 'Cellular processes and response to stress',
    'response to stress'   = 'Cellular processes and response to stress',
    'translation'          = 'Translation and RNA processing',
    'RNA processing'       = 'Translation and RNA processing',
    'transport'            = 'Transport',
    'signal transduction'  = 'Signal transduction and signaling',
    'signaling'            = 'Signal transduction and signaling',
    'proteolysis'          = 'Proteolysis and protein folding',
    'protein folding'      = 'Proteolysis and protein folding'
)

CC <- list(
    'membrane' = 'Membrane',
    'nucleus' = 'Nucleus',
        # 'nuclear protein-containing complex' = 'nuclear protein-containing complex',
        # 'kinetochore' = 'kinetochore',
        # 'nucleosome' = 'nucleosome',
    'cytosol' = 'Cytosol',
    'cytoskeleton' = 'Cytoskeleton',
        # 'spindle' = 'spindle',
    'ribosome' = 'Ribosome',
        # 'translation preinitiation complex' = 'translation preinitiation complex',
    'mitochondrion' = 'Mitochondrion',
    'endomembrane system' = 'Endomembrane system'
        # 'endoplasmic reticulum membrane' = 'endoplasmic reticulum membrane',
        # 'endosome membrane' = 'endosome membrane',
        # 'Golgi membrane' = 'Golgi membrane'
)

MF <- list(
    'catalytic activity'                 = 'Catalytic activity',
    'phosphatase regulator activity'     = 'Catalytic activity',
    'transporter activity'               = 'Transporter activity',
    'structural molecule activity'       = 'Structural molecule activity',
    'binding'                            = 'Binding'
)

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

.getParser <- function() {
    p <- ArgumentParser()
    p$add_argument("-i", help="DE results file (tsv)", type="character", required=TRUE)
    p$add_argument("-g", help="GO results file (ctsv)", type="character", required=TRUE)
    p$add_argument("-o", help="Output directory", type="character", required=TRUE)
    return (p)
}

.isValidArgs <- function(p, args) {
    if (is.null(args$i) | is.null(args$g) | is.null(args$o)) {
        p$print_help()
        stop("Must supply input files and output directory.", call.=FALSE)
    }
}

readDE <- function(f) {
    df <- read.table(f, sep='\t', header=TRUE, quote="")
    df$UniprotAcc <- as.character(df$UniprotAcc)
    df$text <- unlist(lapply(strsplit(df$id, ' \\| '), \(x) x[2]))
    df$id   <- unlist(lapply(strsplit(df$id, ' \\| '), \(x) x[1]))
    df <- df[grep('CNAG', df$id), ]
    df$Significant  <- (abs(df$log2FoldChange) > 1 & df$pvalue < 0.05)
    df <- df[(df$Significant), ]
    return (df)
}

readGO <- function(f, pfx) {
    df <- read.table(f, sep='\t', header=TRUE, quote="")
    return (df)
}

getGOTerms <- function(df, g) {
    df <- df[df$Term %in% names(g), ]
    df <- separate_rows(df, IDs, sep = ",")
    df$Parent <- unlist(g[df$Term])
    return (df)
}

#------------------- Main -----------------------------------#

# Parse command-line arguments
p    <- .getParser()
args <- p$parse_args()
.isValidArgs(p, args)
print(args)

de_df <- readDE(args$i)
de_df <- de_df[, c('id', 'UniprotAcc')]

go_df <- readGO(args$g)
go_df <- go_df[, c('Category', 'GO.ID', 'Term', 'IDs')]
go_df <- getGOTerms(go_df, CC)

res_df <- merge(go_df, de_df, by.x=c('IDs'), by.y=c('UniprotAcc'), all.y=TRUE)
res_df <- res_df[, c('IDs', 'id', 'Category', 'GO.ID', 'Term', 'Parent')]
res_df[is.na(res_df)] <- 'Unknown / hypothetical'
# res_df <- res_df[!is.na(res_df$GO.ID), ]
print(dim(res_df))
print(length(unique(res_df$id)))

## Pie chart
res_df <- res_df %>% group_by(Parent) %>% summarise(count=length(unique(id)), .groups='keep')
print(res_df)
print(res_df %>% summarise(count=sum(count)))
p <- (
    ggplot(res_df, aes(x='', y=count, fill=Parent))
    + geom_bar(stat='identity', position='fill', size=0.1, color='white')
    + scale_fill_brewer(palette = "BrBG")
    # + geom_text(aes(label=x), vjust=0)
    + coord_polar('y', start=0)
    + theme_void()
    + theme(
        legend.title = element_blank(),
        legend.position="bottom"
    )
)
dir.create(args$o, showWarnings=FALSE, recursive=TRUE)
ggsave(file.path(args$o, 'PIE.pdf'), p)
