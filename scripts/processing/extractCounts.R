## Load required libraries
library(mariner)
library(SummarizedExperiment)

## List hicFiles
hicFiles <- list.files(path = "data/raw/hic/replicates", full.names = TRUE)

## Load merged loops
loops <- readRDS("data/mergedLoops.rds")

## Extract replicate Hi-C counts for each .hic file
loopCounts <- pullHicPixels(
    x = loops,
    binSize = 10e3,
    files = hicFiles,
    h5File = "data/mergedLoopCounts.h5",
    norm = "NONE"
)

## Add md5sums (optional)
colData(loopCounts)$md5 <- tools::md5sum(colData(loopCounts)$files)

## Save results
saveRDS(loopCounts, file = "data/mergedLoopCounts.rds")
