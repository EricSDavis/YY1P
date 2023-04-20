## Load required libraries
library(mariner)
library(SummarizedExperiment)
library(DESeq2)

## Load loop counts
loopCounts <- readRDS("data/mergedLoopCounts.rds")

## Update path to count data
path(assay(loopCounts)) <- "data/mergedLoopCounts.h5"

## Build colData from filenames
coldata <- colData(loopCounts)$fileNames |>
    strsplit(split = "_") |>
    do.call(rbind, args = _) |>
    {
        \(x) x[, 3:5]
    }() |> # subset
    `colnames<-`(c("genotype", "biorep", "techrep")) |>
    as.data.frame()
coldata$replicate <- factor(rep(1:4, 2))
coldata$genotype <- factor(coldata$genotype)

## Build DESeq dataset & run DESeq analysis
dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = mariner::counts(loopCounts),
    colData = coldata,
    design = ~ replicate + genotype
)
dds <- DESeq(dds)

## Get shrunken results
res <- lfcShrink(dds, coef = "genotype_WT_vs_KO", type = "apeglm")

## Add to rowData of InteractionMatrix
rowData(loopCounts) <- cbind(rowData(loopCounts), res)

which(res$padj <= 0.05 & res$log2FoldChange > 0) |> length()
which(res$padj <= 0.05 & res$log2FoldChange < 0) |> length()

## Save object with differential results
saveRDS(loopCounts, file = "data/diffLoopCounts.rds")
