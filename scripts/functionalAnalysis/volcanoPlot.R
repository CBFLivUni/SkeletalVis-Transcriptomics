#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))

option_list <- list()
  
option_list$differentialExpression <- make_option('--differentialExpression', type='character',help="path to differential expression table")
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical', help="does the experiment have replicates?")
option_list$volcanoOutput <- make_option('--volcanoOutput', type='character', help="path for output volcano plot")
  
  
opt <- parse_args(OptionParser(option_list=option_list))

#code modified from the enhanced volcano bioconductor package
EnhancedVolcano <- function (toptable, lab, x, y, selectLab = NULL, xlim = c(min(toptable[, x], na.rm = TRUE),
                            max(toptable[, x], na.rm = TRUE)), ylim = c(0, max(-log10(toptable[, y]), na.rm = TRUE) + 5), xlab = bquote(~Log[2] ~"fold change"),
                            ylab = bquote(~-Log[10] ~ italic(P)), axisLabSize = 16, pCutoff = 0.05, pLabellingCutoff = pCutoff, FCcutoff = 2, 
                             title = "", titleLabSize = 16, transcriptPointSize = 0.8, transcriptLabSize = 3, col = c("grey30", "forestgreen", "royalblue","red2"), colOverride = NULL, colAlpha = 1/2, legend = c("NS","Log2 FC", "P", "P & Log2 FC"), legendPosition = "top", 
                             legendLabSize = 10, legendIconSize = 3, DrawConnectors = FALSE, 
                             widthConnectors = 0.5, colConnectors = "black", cutoffLineType = "longdash", 
                             cutoffLineCol = "black", cutoffLineWidth = 0.4, gridlines.major = TRUE, 
                             gridlines.minor = TRUE, border = "partial", borderWidth = 0.8, 
                             borderColour = "black") 
{
  if (!requireNamespace("ggplot2")) {
    stop("Please install ggplot2 first.", call. = FALSE)
  }
  if (!requireNamespace("ggrepel")) {
    stop("Please install ggrepel first.", call. = FALSE)
  }
  if (!is.numeric(toptable[, x])) {
    stop(paste(x, " is not numeric!", sep = ""))
  }
  if (!is.numeric(toptable[, y])) {
    stop(paste(y, " is not numeric!", sep = ""))
  }
  i <- xvals <- yvals <- Sig <- NULL
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(abs(toptable[, x]) > FCcutoff)] <- "FC"
  toptable$Sig[(toptable[, y] < pCutoff)] <- "P"
  toptable$Sig[(toptable[, y] < pCutoff) & (abs(toptable[, 
                                                         x]) > FCcutoff)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig, levels = c("NS", "FC", 
                                                  "P", "FC_P"))
  if (min(toptable[, y], na.rm = TRUE) == 0) {
    warning(paste("One or more P values is 0.", "Converting to minimum possible value..."), 
            call. = FALSE)
    toptable[which(toptable[, y] == 0), y] <- .Machine$double.xmin
  }
  toptable$lab <- lab
  toptable$xvals <- toptable[, x]
  toptable$yvals <- toptable[, y]
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(), 
                                         plot.title = element_text(angle = 0, size = titleLabSize, 
                                                                   face = "bold", vjust = 1), axis.text.x = element_text(angle = 0, 
                                                                                                                         size = axisLabSize, vjust = 1), axis.text.y = element_text(angle = 0, 
                                                                                                                                                                                    size = axisLabSize, vjust = 1), axis.title = element_text(size = axisLabSize), 
                                         legend.position = legendPosition, legend.key = element_blank(), 
                                         legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = legendLabSize), 
                                         title = element_text(size = legendLabSize), legend.title = element_blank())
  if (!is.null(colOverride)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colOverride))), 
                 alpha = colAlpha, size = transcriptPointSize) + 
      scale_color_manual(values = colOverride)
  }
  else {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig)), alpha = colAlpha, 
                 size = transcriptPointSize) + scale_color_manual(values = c(NS = col[1], 
                                                                             FC = col[2], P = col[3], FC_P = col[4]), labels = c(NS = legend[1], 
                                                                                                                                 FC = paste(legend[2], sep = ""), P = paste(legend[3], 
                                                                                                                                                                            sep = ""), FC_P = paste(legend[4], sep = "")))
  }
  plot <- plot + xlab(xlab) + ylab(ylab) + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) + ggtitle(title) + geom_vline(xintercept = c(-FCcutoff, 
                                                                        FCcutoff), linetype = cutoffLineType, colour = cutoffLineCol, 
                                                         size = cutoffLineWidth) + geom_hline(yintercept = -log10(pCutoff), 
                                                                                              linetype = cutoffLineType, colour = cutoffLineCol, size = cutoffLineWidth)
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = borderColour, 
                                                     fill = NA, size = borderWidth))
  }
  else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(size = borderWidth, 
                                                  colour = borderColour), panel.border = element_blank(), 
                         panel.background = element_blank())
  }
  else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  if (DrawConnectors == TRUE && is.null(selectLab)) {
    plot <- plot + geom_text_repel(data = subset(toptable, 
                                                 toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                                 x]) > FCcutoff), aes(label = subset(toptable, 
                                                                                                                                     toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                                                                                                                     x]) > FCcutoff)[, "lab"]), size = transcriptLabSize, 
                                   segment.color = colConnectors, segment.size = widthConnectors, 
                                   vjust = 1.5)
  }
  else if (DrawConnectors == TRUE && !is.null(selectLab)) {
    plot <- plot + geom_text_repel(data = subset(toptable, 
                                                 !is.na(toptable[, "lab"])), aes(label = subset(toptable, 
                                                                                                !is.na(toptable[, "lab"]))[, "lab"]), size = transcriptLabSize, 
                                   segment.color = colConnectors, segment.size = widthConnectors, 
                                   vjust = 1.5)
  }
  else if (DrawConnectors == FALSE && !is.null(selectLab)) {
    plot <- plot + geom_text(data = subset(toptable, !is.na(toptable[, 
                                                                     "lab"])), aes(label = subset(toptable, !is.na(toptable[, 
                                                                                                                            "lab"]))[, "lab"]), size = transcriptLabSize, check_overlap = TRUE, 
                             vjust = 1.5)
  }
  else if (DrawConnectors == FALSE && is.null(selectLab)) {
    plot <- plot + geom_text(data = subset(toptable, toptable[, 
                                                              y] < pLabellingCutoff & abs(toptable[, x]) > FCcutoff), 
                             aes(label = subset(toptable, toptable[, y] < pLabellingCutoff & 
                                                  abs(toptable[, x]) > FCcutoff)[, "lab"]), size = transcriptLabSize, 
                             check_overlap = TRUE, vjust = 1.5)
  }
  return(plot)
}

getVolcanoPlot <- function(data) {
  data <- as.data.frame(data)
  volPlot <- EnhancedVolcano(data,x = "log2FC",
                             y = "padj",
                             lab=data$GeneSymbol,
                             FCcutoff=log2(1.5),
                             pCutoff=0.05)
  
  
  return(volPlot)
}
  
  
  volcanoPlot <- function (differentialExpression, foldChangeOnly=TRUE,volcanoOutput = "volcanoPlot.png") {

 
    differentialExpression <- na.omit(read.delim(differentialExpression))

    colnames(differentialExpression)[1:2] <- c("GeneSymbol","log2FC")
    if (foldChangeOnly == TRUE) {
      return(NULL)
    }  else {
      colnames(differentialExpression)[3] <- c("padj")
      differentialExpression <- differentialExpression %>% 
        group_by(GeneSymbol) %>% slice(which.min(padj)) %>% 
        as.data.frame %>% na.omit()
      }
   
    g <- getVolcanoPlot(differentialExpression)
    ggsave(filename = volcanoOutput, plot = g, device = "png")

  }
  
  opt <- opt[names(opt) != "help"]
  do.call(volcanoPlot, opt)
  

