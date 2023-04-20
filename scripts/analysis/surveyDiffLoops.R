## Load required packages
library(SummarizedExperiment)
library(plotgardener)
library(InteractionSet, include.only = "interactions")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## Source utility functions
source("scripts/utils/customMultiPlot.R")

## Load loops with differential info &
## update path to count data
loops <- readRDS("data/diffLoopCounts.rds")
path(assay(loops)) <- "data/mergedLoopCounts.h5"

## Identify indices differential loops
indices <- which(rowData(loops)$padj <= 0.05 &
    abs(rowData(loops)$log2FoldChange) > 0)

## Define paths to signal track files
atacFiles <- list.files(
    path = "data/raw/signal",
    pattern = "ATAC.*.bw",
    full.names = TRUE
)

chipFiles <- list.files(
    path = "data/raw/signal",
    pattern = "Chip.*.bw",
    full.names = TRUE
)

rnaFiles <- list.files(
    path = "data/raw/signal",
    pattern = "RNA.*.bw",
    full.names = TRUE
)



## Visualization ------------------------------------

## Ensure all interactions are intrachromosomal
stopifnot(all(seqnames1(loops[indices]) ==
    seqnames2(loops[indices])))

pdf(file = "plots/surveyDiffLoops.pdf", width = 3.75, height = 7.25)
for (i in seq_along(indices)) {
    ## Subset for loop
    loop <- loops[indices][i]

    ## Define shared params
    p <- pgParams(

        ## Hi-C params
        assembly = "hg38",
        chrom = seqnames1(loop),
        chromstart = start1(loop) - 100e3,
        chromend = end2(loop) + 100e3,
        resolution = 10e3,
        zrange = c(0, 80),
        norm = "SCALE",

        ## Positional parameters
        x = 0.5,
        y = 0.5,
        width = 3,
        height = 3,
        space = 0.1,

        ## Gene parameters
        gh = 0.5
    )

    ## Create page
    pageCreate(width = 4, height = 7.5, showGuides = FALSE)

    ## Plot Hi-C
    upper <-
        plotHicSquare(
            params = p,
            data = "data/raw/hic/genotype/YY1P_22RV1_WT_inter_30.hic",
            half = "top"
        )
    lower <-
        plotHicSquare(
            params = p,
            data = "data/raw/hic/genotype/YY1P_22RV1_KO_inter_30.hic",
            half = "bottom"
        )
    annoHeatmapLegend(
        plot = upper,
        x = p$x + p$width + p$space,
        y = p$y,
        width = p$space,
        height = p$height * 0.50,
        fontcolor = "black"
    )
    annoPixels(
        plot = upper,
        data = interactions(loop)
    )
    plotText(
        label = "WT",
        fontsize = 10,
        x = p$x + p$space / 2,
        y = p$y + p$space / 2,
        just = c("left", "top")
    )
    plotText(
        label = "KO",
        fontsize = 10,
        x = p$x + p$width - p$space / 2,
        y = p$y + p$height - p$space / 2,
        just = c("right", "bottom")
    )

    ## Plot signal tracks
    ## (ATAC)
    customMultiPlot(
        fn = atacFiles,
        p = p,
        y = "0.1b",
        labs = gsub(".*seq_(.*).bw", "ATAC \\1", atacFiles),
        cols = "forestgreen"
    )

    ## (ChIP - separate scales for H3K27ac & YY1)
    customMultiPlot(
        fn = chipFiles[1:2],
        p = p,
        y = "0.2b",
        labs = gsub(".*seq_(.*).bw", "ChIP \\1", chipFiles[1:2]),
        cols = "steelblue"
    )
    customMultiPlot(
        fn = chipFiles[3],
        p = p,
        y = "0.1b",
        labs = gsub(".*seq_(.*).bw", "ChIP \\1", chipFiles[3]),
        cols = "steelblue"
    )

    # (RNA)
    customMultiPlot(
        fn = rnaFiles,
        p = p,
        y = "0.2b",
        labs = gsub(".*seq_(.*).bw", "RNA \\1", rnaFiles),
        cols = "orange"
    )

    ## Plot genes
    plotGenes(
        params = p,
        y = "0.2b",
        height = p$gh
    )

    ## Genome label
    annoGenomeLabel(
        params = p,
        plot = upper,
        y = "0.2b",
    )
}
dev.off()
