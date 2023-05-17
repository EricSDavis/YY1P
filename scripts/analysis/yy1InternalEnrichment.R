## Quantify internal YY1 signal at 
## gained/lost/static loops to see if
## internal YY1 is preventing larger loop
## formation.

## Load required libraries
library(SummarizedExperiment)
library(InteractionSet)
library(mariner)
library(rtracklayer)
library(ggplot2)
library(patchwork)

## Load table of differential thresholds
tbl <- read.table("tables/diffLoopThresh.txt", header=TRUE)

## Load differential looping results
loopCounts <- readRDS("data/diffLoopCounts.rds")
loops <- interactions(loopCounts)
res <- rowData(loopCounts)

## Path to YY1 bigWig file
yy1File <- "data/raw/signal/Chip_seq_YY1.bw"

## Create plots for each significance threshold
figures <- vector("list", nrow(tbl))
for (i in seq_len(nrow(tbl))) {
  ## Assign loops to each category (static, up, down)  
  st <- loops[which(res$padj > tbl[i,1])]
  up <- loops[which(res$padj <= tbl[i,1] &
                      res$log2FoldChange > tbl[i,2])]
  dn <- loops[which(res$padj <= tbl[i,1] &
                      res$log2FoldChange < -tbl[i,2])]
  
  ## Swap anchors to ensure start1 < start2
  ## and create GRanges of the intervening
  ## ranges (i.e. end1 to start2).
  reg <- lapply(list(st, up, dn), \(x) {
    x <- swapAnchors(x)
    GRanges(seqnames=seqnames1(x),
            ranges=IRanges(start=end1(x),
                           end=start2(x)))
  })
  names(reg) <- c("Static", "Gained", "Lost")
  
  ## Extract YY1 scores from each region
  ## and normalize to the width in
  ## kilobases. Do this in chunks of 15 since its
  ## a lot of data to read in.
  signalPerKb <- lapply(seq_along(reg), \(i) {
    sp <- split(reg[[i]], f=seq(1, length(reg[[i]]), length.out=15))
    scores <- lapply(sp, \(x) {
      import.bw(yy1File, which=x, as="NumericList") |>
        vapply(sum, numeric(1L))
    }) |> do.call(c, args=_) |> unname()
    scores/(width(reg[[i]])/1000)
  })
  n <- vapply(signalPerKb, length, integer(1L))
  names(signalPerKb) <- paste0(c("Static", "Gained", "Lost"), " (n=", n, ")")
  
  ## Visualize YY1 signal per kilobase
  ## among loop classes
  df <- stack(signalPerKb)
  ylim1 <- boxplot.stats(df$values)$stats[c(1, 5)]
  signalBoxPlot <- 
    ggplot(data=df, aes(y=values, x=ind, fill=ind)) +
    geom_boxplot(outlier.colour=NA) +
    ggtitle("YY1 signal inside loops", 
            sprintf("pval <= %s, logFC threshold of %s",
                    tbl[i,1], tbl[i,2])) +
    ylab("YY1 signal per kilobase") +
    xlab("Differential loops (YY1 KO vs. Control)") +
    coord_cartesian(ylim = ylim1*1.05) +
    theme_bw() +
    theme(axis.title.x=element_blank())
  
  signalBoxPlot
  
  ## Visualze loop size distributions
  ## (using gap distances)
  df <- lapply(reg, width) |> stack()
  loopDistributionPlot <- 
    ggplot(data=df, aes(x=log2(values), color=ind)) +
    stat_density(geom = 'line', position = 'identity', na.rm = TRUE) +
    ggtitle("Loop size distributions (gap distances)",
            sprintf("pval <= %s, logFC threshold of %s",
                    tbl[i,1], tbl[i,2])) +
    theme_bw()
  
  ## Include both plots with patchwork
  figure <- signalBoxPlot + loopDistributionPlot
  figures[[i]] <- figure
}

## Write out to pdf
pdf("plots/yy1InternalEnrichment.pdf", width=12, height=4, onefile=TRUE)
print(figures)
dev.off()