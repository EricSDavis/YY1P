## Export files for making deeptools heatmaps

## Load required packages
library(SummarizedExperiment)
library(InteractionSet)
library(plyranges)

## Load differential looping results
loopCounts <- readRDS("data/diffLoopCounts.rds")
loops <- interactions(loopCounts)
res <- rowData(loopCounts)

## Set p-value and lfc thresholds
pthresh <- 0.01
lfcthresh <- 0

## Assign loops to each category (static, up, down)  
st <- loops[which(res$padj > pthresh)]
up <- loops[which(res$padj <= pthresh &
                    res$log2FoldChange > lfcthresh)]
dn <- loops[which(res$padj <= pthresh &
                    res$log2FoldChange < -lfcthresh)]

## Get the unique regions for each set of loops
reg <- lapply(list(st, up, dn), \(x) {
  regions(reduceRegions(x))
})

## Create output filenames
files <- 
  paste0("data/diffLoopBedRegions/",c("static", "gained", "lost"),
       "_p", pthresh, "_lfc" , lfcthresh, ".bed")

## Write out files
mapply(write_bed, x=reg, file=files)
