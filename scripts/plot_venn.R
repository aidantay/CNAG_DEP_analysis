#------------------- Description & Notes --------------------#

# Rscript plot_venn.R -i DE_results.csv -o out

#------------------- Dependencies ---------------------------#

library(argparse)
library(ggvenn)
library(ggplot2)

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

.getParser <- function() {
    p <- ArgumentParser()
    p$add_argument("-i", help="DE results file (csv)", type="character", required=TRUE)
    p$add_argument("-o", help="Output directory", type="character", required=TRUE)
    return (p)
}

.isValidArgs <- function(p, args) {
    if (is.null(args$i) | is.null(args$o)) {
        p$print_help()
        stop("Must supply input file and output directory.", call.=FALSE)
    }
}

#------------------- Main -----------------------------------#

# Parse command-line arguments
p    <- .getParser()
args <- p$parse_args()
.isValidArgs(p, args)
print(args)

## Read data
df <- read.table(args$i, sep=',', header=TRUE)
df$text <- lapply(strsplit(df$Protein.ID, ' \\| '), \(x) x[2])
df <- df[grep('CNAG', df$Protein.ID), ]
df$CBD_vs_Control_Significant   <- (abs(df$CBD_vs_Control_log2.fold.change) > 1 & df$CBD_vs_Control_p.val < 0.05)
df$CBD_V_vs_Control_Significant <- (abs(df$CBD_V_vs_Control_log2.fold.change) > 1 & df$CBD_V_vs_Control_p.val < 0.05)
df$CBD_vs_CBD_V_Significant     <- (abs(df$CBD_vs_CBD_V_log2.fold.change) > 1 & df$CBD_vs_CBD_V_p.val < 0.05)

x <- df[df$CBD_vs_Control_Significant, ]
y <- df[df$CBD_V_vs_Control_Significant, ]
z <- df[df$CBD_vs_CBD_V_Significant, ]
l <- list(
    'CBD'=x$text,
    'CBDV'=y$text,
    'CBD_CBDV'=z$text
)

p <- ggvenn(l, c('CBD', 'CBDV'), fill_color=c('#1f78b4', '#e31a1c'), text_size=5, auto_scale=TRUE)
dir.create(args$o, showWarnings=FALSE, recursive=TRUE)
ggsave(file.path(args$o, 'VENN.pdf'), p)