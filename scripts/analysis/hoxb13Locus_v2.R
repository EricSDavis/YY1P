library(plotgardener)
library(mariner)
library(SummarizedExperiment)
library(InteractionSet)
library(grid)

## Source utility functions
source("scripts/utils/normalizeHic.R")
source("scripts/utils/customMultiPlot.R")

## Input path of hicFiles
hicFiles <- c(
  "WT" = "data/raw/hic/genotype/YY1P_22RV1_WT_inter_30.hic",
  "KO" = "data/raw/hic/genotype/YY1P_22RV1_KO_inter_30.hic"
)

## Load loops with differential info &
## update path to count data
loops <- readRDS("data/diffLoopCounts.rds")
path(assay(loops)) <- "data/mergedLoopCounts.h5"


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

ctcfFiles <- list.files(
  path = "data/raw/signal",
  pattern = "CTCF.*.bigWig",
  full.names = TRUE
)

## Visualization -----------------------------------------------

## Define region of interest
p <- pgParams(
  
  ## Hi-C params
  assembly = "hg38",
  chrom = "chr17",
  chromstart = 48702638-110e3,
  chromend = 48852638+110e3,
  resolution = 5e3,
  zrange = c(0, 100),
  norm = "SCALE",
  
  ## Positional parameters
  x = 0.5,
  y = 0.5,
  width = 5,
  height = 1.5,
  space = 0.1,
  
  ## Gene parameters
  gh = 0.5,
  
  ## Loop arch parameters
  archHeight = 0.5
)

## Create page
pdf("plots/hoxb13Locus_v2.pdf", width=6, height=8.5)
pageCreate(width=6, height=8.5, showGuides=FALSE)

## Define additional data to pull
## for rectangular plot (logic from
## plotgardener::plotHicRectangle)
dimRatio <- p$height / p$width
buffer <- ((p$chromend - p$chromstart)*dimRatio + p$resolution)

## Normalize Hi-C for differences
## in read-depth (use buffer for
## pulling the right amount of data)
hicData <- normalizeHic(
  files=hicFiles,
  params=p,
  chromstart=p$chromstart-buffer,
  chromend=p$chromend+buffer
)

## Plot Hi-C
wtHic <- plotHicRectangle(
  params=p,
  data=hicData[['WT']]
)
koHic <- plotHicRectangle(
  params=p,
  data=hicData[['KO']],
  y=paste0(p$space, 'b')
)
lapply(list(wtHic, koHic), \(hic) {
  annoHeatmapLegend(
    plot = hic,
    x = p$x + p$width + p$space,
    y = hic$y,
    width = p$space,
    height = p$height * 0.75,
    fontcolor = "black"
  )
})
inset <- unit(p$space/2, "inches")
plotText(
  label=names(hicData),
  x=c(wtHic$x+inset, koHic$x+inset),
  y=c(wtHic$y+inset, koHic$y+inset),
  just=c('left', 'top')
)

## Plot loops
plotPairsArches(
  data=interactions(loops),
  params=p,
  flip=TRUE,
  height=p$archHeight,
  y='0b'
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

## (ChIP - separate scales for H3K27ac & YY1 & CTCF)
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
customMultiPlot(
  fn = ctcfFiles,
  p = p,
  y = "0.1b",
  labs = gsub(".*signal/(.*)_Ch.*", "ChIP \\1", ctcfFiles),
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
  y = paste0(p$space*2, 'b'),
  height = p$gh
)

## Genome label
annoGenomeLabel(
  params = p,
  plot = wtHic,
  at = seq(p$chromstart, p$chromend, by=10e3),
  y = paste0(p$space*2, 'b')
)
dev.off()