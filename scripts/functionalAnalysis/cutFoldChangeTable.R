#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list()

option_list$foldChangeTable <- make_option('--foldChangeTable', type='character',help = "path to foldchange table")
option_list$comparisonNumber <- make_option('--comparisonNumber', type='numeric',help="comparison number to extract")
option_list$cutTable <- make_option('--cutTable', type='character',help="path to output table")


opt <- parse_args(OptionParser(option_list=option_list,description = "Splits a merged differential expression table to produce a table of gene IDs, fold changes and pvalues"))

getColumns <- function(FCTable, columnIndex, numberOfColumns) {
    cbind(FCTable[1], FCTable[seq(columnIndex, columnIndex + numberOfColumns - 1)])
}

cutFoldChangeTable <- function (foldChangeTable, comparisonNumber, cutTable = "cutTable.txt") {
    
    
    foldChangeTable = read.table(foldChangeTable, header = TRUE, sep = "\t")
    columnNames = names(foldChangeTable)
    if (length(grep("P.Val|padj", columnNames)) == 0) {
        returnTable = getColumns(foldChangeTable, comparisonNumber + 1, 1)
    }   else {
        returnTable = getColumns(foldChangeTable, comparisonNumber * 2, 2)
    }
    write.table(returnTable, file = cutTable, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}


opt <- opt[names(opt) != "help"]
do.call(cutFoldChangeTable, opt)


