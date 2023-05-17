## Are gained/lost/static loop anchors
## enriched for CTCF signal?
## As control for the YY1 version

## Load required libraries
library(SummarizedExperiment)
library(InteractionSet)
library(rtracklayer)
library(gridExtra)
library(grid)

## Load table of differential thresholds
tbl <- read.table("tables/diffLoopThresh.txt", header=TRUE)

## Load differential looping results
loopCounts <- readRDS("data/diffLoopCounts.rds")
loops <- interactions(loopCounts)
res <- rowData(loopCounts)

## Path to YY1 bigWig file
ctcfFile <- "data/raw/signal/CTCF_Chipseq_signal_p-value_GRCh38.bigWig"

## Plot YY1 signal for diff loops for each significance
## and logFC threshold
pdf(file="plots/ctcfEnrichment.pdf", width=12, height=4)
for (i in seq_len(nrow(tbl))) {

  ## Assign loops to each category (static, up, down)  
  st <- loops[which(res$padj > tbl[i,1])]
  up <- loops[which(res$padj <= tbl[i,1] &
                      res$log2FoldChange > tbl[i,2])]
  dn <- loops[which(res$padj <= tbl[i,1] &
                      res$log2FoldChange < -tbl[i,2])]
  
  ## Get the unique regions for each set of loops
  reg <- lapply(list(st, up, dn), \(x) {
    regions(reduceRegions(x))
  })
  
  ## Extract scores
  scores <- lapply(seq_along(reg), \(i) {
    import.bw(ctcfFile, which=reg[[i]], as="NumericList") |>
      vapply(sum, numeric(1L))
  })
  n <- vapply(scores, length, integer(1L))
  names(scores) <- paste0(c("Static", "Gained", "Lost"), " (n=", n, ")")
  
  ## Calculate stats
  df <- scores
  names(df) <- c("Static", "Gained", "Lost")
  df <- stack(df) |> setNames(c("vals", "grp"))
  stats <- pairwise.wilcox.test(
    x=df$vals,
    g=df$grp,
    p.adjust.method='bonferroni'
  )
  pvals <- formatC(
    x=stats$p.value,
    format = "e",
    digits = 2
  )
  
  ## Visualize
  title <- sprintf("pval <= %s, logFC threshold of %s",
                   tbl[i,1], tbl[i,2])
  ## Boxplot
  par(mfrow=c(1, 3))
  boxplot(
    x=scores,
    ylab="CTCF signal",
    xlab="Differential Loops (YY1 KO vs. Control)",
    main=title,
    outline=FALSE
  )
  
  ## Violin plot
  vioplot::vioplot(
    x=scores,
    ylab="CTCF signal",
    xlab="Differential Loops (YY1 KO vs. Control)",
    main=title
  )
  
  ## Stats
  vp <- viewport(x=0.75, y=0.5, width=1)
  tg <- tableGrob(pvals, vp=vp, theme=ttheme_minimal())
  grid.draw(tg)
  
  par(mfrow=c(1,1))
  
}
dev.off()



