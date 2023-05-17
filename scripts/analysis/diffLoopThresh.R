## Make a table showing how the number
## of differential loops identified 
## changes with different significance
## and log2FC thresholds.

## Load required libraries
library(SummarizedExperiment)

## Load differential looping data
loopCounts <- readRDS("data/diffLoopCounts.rds")

## Get differential results from rowData
res <- rowData(loopCounts)

## Iterate through different combos of
## padj and log2FC thresholds
pvthresh <- c(0.01, 0.05, 0.1)
fcthresh <- c(0, 1, 1.5, 2)
grid <- expand.grid(pvthresh=pvthresh, fcthresh=fcthresh)

for (i in seq_len(nrow(grid))) {
  grid$up[i] <- which(res$padj <= grid[i,1] &
              res$log2FoldChange > grid[i,2]) |> length()
  grid$dn[i] <- which(res$padj <= grid[i,1] &
              res$log2FoldChange < -grid[i,2]) |> length()
  grid$total[i] <- grid$up[i] + grid$dn[i]
}

## Write out results
write.table(grid, file="tables/diffLoopThresh.txt", sep="\t",
            row.names=FALSE)
