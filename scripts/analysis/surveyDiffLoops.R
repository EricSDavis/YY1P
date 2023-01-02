## Load required packages
library(SummarizedExperiment)
library(plotgardener)
library(InteractionSet, include.only = "interactions")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## Load loops with differential info &
## update path to count data
loops <- readRDS("data/diffLoopCounts.rds")
path(assay(loops)) <- "data/mergedLoopCounts.h5"

## Identify indices differential loops
indices <- which(rowData(loops)$padj <= 0.05 &
    abs(rowData(loops)$log2FoldChange) > 0)



## Visualization ------------------------------------

## Ensure all interactions are intrachromosomal
stopifnot(all(seqnames1(loops[indices]) ==
    seqnames2(loops[indices])))

pdf(file = "plots/surveyDiffLoops.pdf", width = 4, height = 5)
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
        norm = "KR",

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
    pageCreate(width = 4, height = 5, showGuides = FALSE)

    ## Plot Hi-C
    upper <-
        plotHicSquare(
            params = p,
            data = "data/raw/hic/condition/WT_inter_30.hic",
            half = "top"
        )
    lower <-
        plotHicSquare(
            params = p,
            data = "data/raw/hic/condition/KO_inter_30.hic",
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

    ## Plot genes
    plotGenes(
        params = p,
        y = p$y + p$height + p$space,
        height = p$gh
    )

    ## Genome label
    annoGenomeLabel(
        params = p,
        plot = upper,
        y = p$y + p$height + p$gh + p$space,
    )
}
dev.off()
