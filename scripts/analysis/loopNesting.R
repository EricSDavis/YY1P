## Testing loop nesting.
## Are lost loops (dn) inside gained loops (up)
## after KO of YY1?

## Load required libraries
library(InteractionSet)
library(SummarizedExperiment)
library(ggplot2)
library(patchwork)

## Load table of differential thresholds
tbl <- read.table("tables/diffLoopThresh.txt", header=TRUE)

## Load differential looping results
loopCounts <- readRDS("data/diffLoopCounts.rds")
loops <- interactions(loopCounts)
res <- rowData(loopCounts)

## Create plots for each significance threshold
figures <- vector("list", nrow(tbl))
for (i in seq_len(nrow(tbl))) {
  ## Assign loops to each category (all, static, up, down)  
  all <- loops
  st <- loops[which(res$padj > tbl[i,1])] #not used
  up <- loops[which(res$padj <= tbl[i,1] &
                      res$log2FoldChange > tbl[i,2])]
  dn <- loops[which(res$padj <= tbl[i,1] &
                      res$log2FoldChange < -tbl[i,2])]
  sets <- list("All"=all, "Gained"=up, "Lost"=dn)
  
  ## Compute the inner arc of each loop
  inner <- lapply(sets, \(x) {
    x <- swapAnchors(x) # ensure start1 < start2
    GRanges(seqnames=seqnames1(x),
            ranges=IRanges(start=end1(x),
                           end=start2(x)))
  })
  
  ## Calculate overlap functions (x within y)
  ## and as a percentage of total x
  calcOverlap <- \(x,y) {
    sum(x %within% y)
  }
  calcOverlapPerc <- \(x,y) {
    calcOverlap(x,y)/length(x) * 100
  }
  
  ## Compute percent of x within y (rows within cols)
  outer(inner, inner, FUN=\(x, y) {
    mapply(calcOverlapPerc, x=x, y=y) |>
      round(2)
  })
  
  ## Numbers (not percent)
  outer(inner, inner, FUN=\(x, y) {
    mapply(calcOverlap, x=x, y=y) |>
      round(2)
  })
  
  ## How often do we see a gained loop inside
  ## any loop by chance?
  set.seed(123)
  gil <- vapply(seq_len(1000), \(x) {
    s <- sample(x=length(inner$All),
                size=length(inner$Gained),
                replace=TRUE)
    calcOverlap(inner$Gained, inner$All[s])
  }, numeric(1L))
  
  gdf <- data.frame(value=gil)
  gilPlot <- ggplot(data=gdf, aes(x=value)) +
    stat_density(geom="line") +
    labs(title="Gained loops inside lost loops",
         subtitle=sprintf("pval <= %s, logFC threshold of %s",
                          tbl[i,1], tbl[i,2]),
         x="Number of overlaps",
         y="Density") +
    geom_vline(xintercept=calcOverlap(inner$Gained, inner$Lost),
               color="blue") +
    geom_text(aes(x=calcOverlap(inner$Gained, inner$Lost),
                  y=max(density(gdf$value)$y),
                  label="Observed overlap"),
              hjust=-0.1,
              color="blue") +
    theme_bw()
  
  ## How often do we see a lost loop inside
  ## any loop by chance?
  set.seed(123)
  lig <- vapply(seq_len(1000), \(x) {
    s <- sample(x=length(inner$All),
                size=length(inner$Lost),
                replace=TRUE)
    calcOverlap(inner$Lost, inner$All[s])
  }, numeric(1L))
  
  ldf <- data.frame(value=lig)
  ligPlot <- ggplot(data=ldf, aes(x=value)) +
    stat_density(geom="line") +
    labs(title="Lost loops inside gained loops",
         subtitle=sprintf("pval <= %s, logFC threshold of %s",
                          tbl[i,1], tbl[i,2]),
         x="Number of overlaps",
         y="Density") +
    geom_vline(xintercept=calcOverlap(inner$Lost, inner$Gained),
               color="blue") +
    geom_text(aes(x=calcOverlap(inner$Lost, inner$Gained),
                  y=max(density(ldf$value)$y),
                  label="Observed overlap"),
              hjust=-0.1,
              color="blue") +
    theme_bw()
  
  ## Include both plots with patchwork
  figure <- gilPlot + ligPlot
  figures[[i]] <- figure
}

## Write out to pdf
pdf("plots/loopNesting.pdf", width=12, height=4, onefile=TRUE)
print(figures)
dev.off()

