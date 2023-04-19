#########################################
## Merge SIP loops called in bioreps & ##
## Condition-merged datasets           ##
#########################################

## Load required packages
library(mariner)
library(data.table)

## List loop files
loopFiles <- list.files(path="data/raw/loops",
                        full.names=TRUE)

## Read in loop files & convert to GInteractions list
loopList <- 
  lapply(loopFiles, fread) |>
  lapply(as_ginteractions) |>
  setNames(gsub("^.*loops/(.*)_inter.*$", "\\1", loopFiles))

## Cluster, merge & bin at 10e3
loops <- 
  mergePairs(x=loopList,
             radius=10e3,
             column="APScoreAvg") |>
  binPairs(binSize=10e3)

## Save output
saveRDS(object=loops, file="data/mergedLoops.rds")
