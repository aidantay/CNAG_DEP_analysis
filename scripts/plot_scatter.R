#------------------- Description & Notes --------------------#

# Rscript plot_scatter.R -i DE_results.csv -o out

#------------------- Dependencies ---------------------------#

library(argparse)
library(ggplot2)
library(ggrepel)
library(tidyverse)

#------------------- Constants ------------------------------

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

.getParser <- function() {
    p <- ArgumentParser()
    p$add_argument("-i", help="DE results file (tsv)", type="character", required=TRUE)
    p$add_argument("-o", help="Output directory", type="character", required=TRUE)
    return (p)
}

.isValidArgs <- function(p, args) {
    if (is.null(args$i) | is.null(args$o)) {
        p$print_help()
        stop("Must supply input files and output directory.", call.=FALSE)
    }
}


plot <- function(df, x, y) {
    x_lfc_col  <- paste(x, 'log2.fold.change', sep='_')
    x_pval_col <- paste(x, 'p.val', sep='_')
    y_lfc_col  <- paste(y, 'log2.fold.change', sep='_')
    y_pval_col <- paste(y, 'p.val', sep='_')

    cond1 = abs(df[[x_lfc_col]]) > 1
    cond2 = abs(df[[y_lfc_col]]) > 1
    cond3 = df[[x_pval_col]] < 0.05
    cond4 = df[[y_pval_col]] < 0.05
    
    df$type <- "Not significant"
    df$type[cond1 & cond3] <- "Significant in CBD or CBDV"
    df$type[cond2 & cond4] <- "Significant in CBD or CBDV"
    df$type[cond1 & cond2 & cond3 & cond4] <- "Significant in CBD and CBDV"

    ids <- c(
        'CNAG_02830', 
        'CNAG_03196',
        'CNAG_01850', 'CNAG_00130',
        'CNAG_05925', 'CNAG_03993', 'CNAG_05172',
        'CNAG_00730', 'CNAG_00869', 'CNAG_00895',
        'CNAG_01216', 'CNAG_04753'
    )    
    df$text[!df$text %in% ids] <- ""

    cols <- c('Significant in CBD and CBDV'='#3f007d',
              'Significant in CBD or CBDV'='#807dba',
              'Not significant'='#e0e0e5'
    )    
    p <- (
        ggplot(data=df, aes(x=.data[[x_lfc_col]], y=.data[[y_lfc_col]], col=type))
        + geom_point()
        + geom_text_repel(aes(label=text), size=5, max.overlaps=500, min.segment.length=0.1)
        + scale_color_manual(values=cols)
        + geom_vline(xintercept=0, linetype='dashed', color='grey')
        + geom_hline(yintercept=0, linetype='dashed', color='grey')
        + labs(
            x=expression('CBD (Log'['2']*'FC)'),
            y=expression('CBDV (Log'['2']*'FC)')
        )
        + theme_classic()
        + theme(
            text = element_text(size = 20),
            legend.title = element_blank(),
            legend.position = c(0.01, 0.99), 
            legend.justification = c(0, 1))
    )
    return(p)
}

#------------------- Main -----------------------------------#

# Parse command-line arguments
p    <- .getParser()
args <- p$parse_args()
.isValidArgs(p, args)
print(args)

df <- read.table(args$i, sep=',', header=TRUE)
df$text <- lapply(strsplit(df$Protein.ID, ' \\| '), \(x) x[2])
df <- df[grep('CNAG', df$Protein.ID), ]
contrasts <- c('CBD_vs_Control', 'CBD_V_vs_Control', 'CBD_vs_CBD_V')
contrasts <- combn(contrasts,2)
plots <- apply(contrasts, 2, \(x) plot(df, x[1], x[2]))

dir.create(args$o, showWarnings=FALSE, recursive=TRUE)
ggsave(file.path(args$o, 'SCATTER.pdf'), plots[[1]], width=10, height=10)