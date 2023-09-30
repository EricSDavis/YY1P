## Test aggreTAD with mariner
## on some real data

## Load required libraries
library(SummarizedExperiment)

## Load differential looping data
loopCounts <- readRDS("data/diffLoopCounts.rds")
res <- rowData(loopCounts)


## Write a helper function to convert
## loops into TADs
makeTadsFromLoops <- function(x) {
  d <- InteractionSet::pairdist(x)/2
  df <- data.frame(
    seqnames1(x),
    start1(x)-d,
    end2(x)+d
  )
  as_ginteractions(cbind(df, df))
}

tads <- makeTadsFromLoops(interactions(loopCounts))

## Separate into gained/lost
gained <- tads[which(res$padj < 0.05 & res$log2FoldChange > 0)]
lost <- tads[which(res$padj < 0.05 & res$log2FoldChange < 0)]

## Path to .hic files
ko <- "data/raw/hic/genotype/YY1P_22RV1_KO_inter_30.hic"
wt <- "data/raw/hic/genotype/YY1P_22RV1_WT_inter_30.hic"

gainedTads <- pullHicMatrices(
  x=gained,
  files=wt,
  binSize=10e3,
  blockSize=100e6
)

reg <- regularize(gainedTads, ndim=c(100, 100), nBlocks=10, scale=TRUE)

mat <- aggHicMatrices(reg)

# plotMatrix(mat, zrange=c(0, 100000))
plotMatrix(mat, zrange=c(0, 0.1))



#############################


library(raster)
library(tidyverse)
library(data.table)
library(reshape2)
library(ggplot2)
options(scipen=999)

########################################################
#################### FUNCTIONS #########################
########################################################

#' @param loops data frame of regions to aggregate. Fist 6 columns must be in BEPDE format
#' @param hic  path to .hic file to use for aggregation
#' @param buffer fraction of the loops to add to each side (e.g. 0.5)
#' @param res resolution of hic data to use
#' @param size integer descibing the number of bins in the aggregated matix. (e.g. a size of 100 would result in an aggregated matrix with dimensions 100 x 100)

### Define a function that builds aggregate TADS
aggregateTAD <- function(loops,hic,buffer=.5,res,size=100)
{
  # remove interchrom loops
  loops = loops[which(loops[,1] == loops[,4]),]
  
  # define window to plot
  loops$size = (loops[,5]- loops[,2])/res
  loops$bufferstart = res*(loops[,2]/res - round(loops$size*buffer))
  loops$bufferend = res*(loops[,5]/res + round(loops$size*buffer))
  loops$coords = paste0(loops[,1],":",loops$bufferstart,":",loops$bufferend)
  
  # create matrix to fill
  aggreTAD = matrix(0,nrow=size,ncol=size)
  
  # keep track of total counts
  totalcounts = 0
  
  # iterate through loops
  for (i in 1:nrow(loops))
  {
    # print update
    print (paste(i, "of",nrow(loops),"loops"))
    
    # get loop info
    loop = loops[i,]
    
    # get pixels
    sparseMat = as.data.table(strawr::straw("KR", hic, loop$coords, loop$coords, "BP", res))
    
    # define bins
    startcoord = loop$bufferstart
    endcoord   = loop$bufferend
    bins <- seq(from = startcoord, to = endcoord, by = res)
    
    # make empty long format matrix
    longmat = as.data.table(expand.grid(bins,bins))
    longmat$counts = 0
    colnames(longmat) = c("x","y","counts")
    
    ## Set keys
    setkeyv(sparseMat, c('x', 'y'))
    
    ## Get counts by key
    longmat$counts <- sparseMat[longmat]$counts
    
    ## Set unmatched counts (NA) to 0
    longmat[is.na(counts), counts := 0]
    
    # convert to wide matrix
    wideMat <- reshape2::acast(longmat, x ~ y, value.var = 'counts')
    
    # make symmetric
    wideMat[lower.tri(wideMat)] = t(wideMat)[lower.tri(wideMat)]
    
    # resize the matrix
    r_wideMat <- raster(wideMat)
    extent(r_wideMat) <- extent(c(-180, 180, -90, 90))
    
    resizedMat <- raster(ncol=size,  nrow=size)
    resizedMat <- resample(r_wideMat, resizedMat)
    
    # convert to matrix
    resizedMat = as.matrix(resizedMat)
    # resizedMat[is.na(resizedMat)] <- 0
    
    # update counts
    currentcounts = sum(resizedMat,na.rm = TRUE)
    totalcounts = totalcounts + currentcounts
    
    # normalize to more evenly weight different sized loops
    resizedMatNorm = resizedMat / currentcounts
    
    aggreTAD = aggreTAD + resizedMatNorm
  }
  
  # normalize for counts and number of loops
  aggreTAD = aggreTAD*totalcounts/nrow(loops)
  
  return (aggreTAD)
}


#' @param AggTAD Aggregated TAD object to plot
#' @param maxval integer representing the maximum value to plot (essentialy sets the top of the zrange)
#' @param cols color palette for plotting
#' @param title title of the plot

### Define a function that builds aggregate TADS
plotAggTAD <- function(AggTAD,maxval = 10000,cols = RColorBrewer::brewer.pal(6,"YlGnBu"),title="")
{
  # Convert to long format for ggplot
  AggTAD_long = setNames(melt(AggTAD), c('x', 'y', 'counts'))
  AggTAD_long$counts = AggTAD_long$counts
  
  ggplot(data=AggTAD_long,mapping=aes(x=x,y=y,fill=counts)) + 
    geom_tile() + 
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(aspect.ratio = 1) + 
    ggtitle(title) +
    scale_fill_gradientn(colours = cols,
                         na.value=cols[maxval],
                         limits=c(0,maxval),
                         oob = scales::squish) 
}


########################################################
####################### CODE ###########################
########################################################

#---------------------- CODE LIMA ----------------------#

library(data.table)
library(raster)
library(ggplot2)
# set hic file
hicfile = ko

loops <- interactions(loopCounts)[which(res$padj < 0.05 & res$log2FoldChange > 0)]
# get loops and sizes (in bins)
loops <- as.data.frame(loops)[,c(1:3,6:8)]

# Aggregate TADs
metaTAD = aggregateTAD(loops=loops[1:50,],
                       hic=hicfile,
                       buffer=0.5,
                       res=10000,
                       size=100)


# Convert to long format for ggplot
metaTAD_long = setNames(melt(metaTAD), c('x', 'y', 'counts'))
metaTAD_long$counts = metaTAD_long$counts

# Plot
cols = RColorBrewer::brewer.pal(6,"YlGnBu")
maxval = 2500
ggplot(data=metaTAD_long,mapping=aes(x=x,y=y,fill=counts)) + 
  geom_tile() + 
  theme_void() + 
  theme(aspect.ratio = 1) + 
  scale_fill_gradientn(colours = cols,
                       na.value=cols[maxval],
                       limits=c(0,maxval),
                       oob = scales::squish) 


#---------------------- CODE YAPP (TAD/Loop loss) ----------------------#

# set hic file
hicfile_control  = "/Users/dphansti/Dropbox (Personal)/Work/Projects/Ongoing/YAP/data/hg19/hic/HEK/YAPP_HEK_control_inter_30.hic"
hicfile_sorbitol = "/Users/dphansti/Dropbox (Personal)/Work/Projects/Ongoing/YAP/data/hg19/hic/HEK/YAPP_HEK_sorbitol_inter_30.hic"

# get loops and sizes (in bins)
loops_control = "/Users/dphansti/Dropbox (Personal)/Work/Projects/Ongoing/YAP/data/hg19/DiffLoopCalls/YAP_control_5kbLoops.txt"
loops = data.frame(read_tsv(loops_control,col_names = TRUE))
loops = loops[which(loops$chromosome1 == loops$chromosome2),]
loops$chromosome1 = gsub("chr","",loops$chromosome1)
loops$chromosome2 = gsub("chr","",loops$chromosome2)

# Aggregate TADs
metaTAD_control = aggregateTAD(loops=loops[1:1000,],
                               hic=hicfile_control,
                               buffer=0.5,
                               res=5000,
                               size=100)

metaTAD_sorbitol = aggregateTAD(loops=loops[1:1000,],
                                hic=hicfile_sorbitol,
                                buffer=0.5,
                                res=5000,
                                size=100)

# Convert to long format for ggplot
metaTAD_control_long = setNames(melt(metaTAD_control), c('x', 'y', 'counts'))
metaTAD_control_long$counts = metaTAD_control_long$counts

metaTAD_sorbitol_long = setNames(melt(metaTAD_sorbitol), c('x', 'y', 'counts'))
metaTAD_sorbitol_long$counts = metaTAD_sorbitol_long$counts

# Plot
cols = RColorBrewer::brewer.pal(6,"YlGnBu")
maxval = 1500
ggplot(data=metaTAD_control_long,mapping=aes(x=x,y=y,fill=counts)) + 
  geom_tile() + 
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(aspect.ratio = 1) + 
  ggtitle("control loops in untreated cells") +
  scale_fill_gradientn(colours = cols,
                       na.value=cols[maxval],
                       limits=c(0,maxval),
                       oob = scales::squish) 

ggplot(data=metaTAD_sorbitol_long,mapping=aes(x=x,y=y,fill=counts)) + 
  geom_tile() + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(aspect.ratio = 1) +
  ggtitle("control loops in treated cells") +
  scale_fill_gradientn(colours = cols,
                       na.value=cols[maxval],
                       limits=c(0,maxval),
                       oob = scales::squish) 
