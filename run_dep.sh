#!/bin/sh

################### Workspace & Notes #########################

################### Dependencies ##############################

set -e

################### Global Variables ##########################

################### Functions #################################

############################ Main #############################

Rscript scripts/enrichment.R \
    -i CBD_V_vs_Control.tsv \
    -u GO_KEGG/uniprot/geneID2GO_all.tsv \
    -o out \
    -p

Rscript scripts/enrichment.R \
    -i CBD_vs_Control.tsv \
    -u GO_KEGG/uniprot/geneID2GO_all.tsv \
    -o out \
    -p

Rscript scripts/plot_venn.R \
    -i DE_results.csv \
    -o out

Rscript scripts/plot_pie.R \
    -i CBD_V_vs_Control.tsv \
    -g GO_KEGG/uniprot/all/GO_All_CBD_V_vs_Control.tsv \
    -o out

Rscript scripts/plot_pie.R \
    -i CBD_vs_Control.tsv \
    -g GO_KEGG/uniprot/all/GO_All_CBD_vs_Control.tsv \
    -o out

Rscript scripts/plot_scatter.R \
    -i DE_results.csv \
    -o out