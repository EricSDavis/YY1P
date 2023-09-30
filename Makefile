.PHONY: clean

objects :=\
	data/mergedLoops.rds\
	data/mergedLoopCounts.h5\
	data/mergedLoopCounts.rds\
	data/diffLoopCounts.rds\
	tables/diffLoopThresh.txt\
	plots/surveyDiffLoops.pdf\
	plots/hoxb13Locus.pdf\
	plots/hoxb13Locus_v2.pdf\
	plots/yy1Enrichment.pdf\
	plots/ctcfEnrichment.pdf\
	data/diffLoopBedRegions/static_p0.01_lfc0.bed\
	data/diffLoopBedRegions/gained_p0.01_lfc0.bed\
	data/diffLoopBedRegions/lost_p0.01_lfc0.bed\
	plots/yy1InternalEnrichment.pdf\
	plots/loopNesting.pdf

all: $(objects)

clean:
	rm -rf $(objects)
	
###################################
## Differential Looping Analysis ##
###################################

## Merge loops
data/mergedLoops.rds:\
	data/raw/loops/YY1P_22RV1_KO_inter_30_5kbLoops.txt\
	data/raw/loops/YY1P_22RV1_WT_inter_30_5kbLoops.txt\
	data/raw/loops/YY1P_22RV1_KO_1_1_inter_30_5kbLoops.txt\
	data/raw/loops/YY1P_22RV1_KO_1_2_inter_30_5kbLoops.txt\
	data/raw/loops/YY1P_22RV1_KO_2_1_inter_30_5kbLoops.txt\
	data/raw/loops/YY1P_22RV1_KO_2_2_inter_30_5kbLoops.txt\
	data/raw/loops/YY1P_22RV1_WT_1_1_inter_30_5kbLoops.txt\
	data/raw/loops/YY1P_22RV1_WT_1_2_inter_30_5kbLoops.txt\
	data/raw/loops/YY1P_22RV1_WT_2_1_inter_30_5kbLoops.txt\
	data/raw/loops/YY1P_22RV1_WT_2_2_inter_30_5kbLoops.txt\
	scripts/processing/mergeLoops.R
		mkdir -p data
		Rscript scripts/processing/mergeLoops.R 

## Extract counts
data/mergedLoopCounts.h5\
data/mergedLoopCounts.rds:\
	data/raw/hic/replicates/YY1P_22RV1_KO_1_1_inter_30.hic\
	data/raw/hic/replicates/YY1P_22RV1_KO_1_2_inter_30.hic\
	data/raw/hic/replicates/YY1P_22RV1_KO_2_1_inter_30.hic\
	data/raw/hic/replicates/YY1P_22RV1_KO_2_2_inter_30.hic\
	data/raw/hic/replicates/YY1P_22RV1_WT_1_1_inter_30.hic\
	data/raw/hic/replicates/YY1P_22RV1_WT_1_2_inter_30.hic\
	data/raw/hic/replicates/YY1P_22RV1_WT_2_1_inter_30.hic\
	data/raw/hic/replicates/YY1P_22RV1_WT_2_2_inter_30.hic\
	data/mergedLoops.rds\
	scripts/processing/extractCounts.R
		mkdir -p data
		rm -rf data/mergedLoopCounts.h5
		Rscript scripts/processing/extractCounts.R

## Differential analysis
data/diffLoopCounts.rds:\
	data/mergedLoopCounts.h5\
	data/mergedLoopCounts.rds\
	scripts/analysis/differentialLoops.R
		mkdir -p data
		Rscript scripts/analysis/differentialLoops.R
		
## Apply significance thresholds
tables/diffLoopThresh.txt:\
	data/diffLoopCounts.rds\
	scripts/analysis/diffLoopThresh.R
		mkdir -p tables
		Rscript scripts/analysis/diffLoopThresh.R

## Survey of top diff loops
plots/surveyDiffLoops.pdf:\
	data/diffLoopCounts.rds\
	data/mergedLoopCounts.h5\
	data/raw/hic/genotype/YY1P_22RV1_WT_inter_30.hic\
	data/raw/hic/genotype/YY1P_22RV1_KO_inter_30.hic\
	data/raw/signal/ATAC_seq_EV.bw\
	data/raw/signal/ATAC_seq_YY1KD.bw\
	data/raw/signal/Chip_seq_H3K27ac_.bw\
	data/raw/signal/Chip_seq_INPUT.bw\
	data/raw/signal/Chip_seq_YY1.bw\
	data/raw/signal/RNA_seq_EV.bw\
	data/raw/signal/RNA_seq_YY1KD.bw\
	scripts/utils/customMultiPlot.R\
	scripts/analysis/surveyDiffLoops.R
		mkdir -p plots
		Rscript scripts/analysis/surveyDiffLoops.R


##############################################
## Make plot of region of interest (HOXB13) ##
##############################################

plots/hoxb13Locus.pdf:\
	data/diffLoopCounts.rds\
	data/mergedLoopCounts.h5\
	data/raw/hic/genotype/YY1P_22RV1_WT_inter_30.hic\
	data/raw/hic/genotype/YY1P_22RV1_KO_inter_30.hic\
	data/raw/signal/ATAC_seq_EV.bw\
	data/raw/signal/ATAC_seq_YY1KD.bw\
	data/raw/signal/Chip_seq_H3K27ac_.bw\
	data/raw/signal/Chip_seq_INPUT.bw\
	data/raw/signal/Chip_seq_YY1.bw\
	data/raw/signal/RNA_seq_EV.bw\
	data/raw/signal/RNA_seq_YY1KD.bw\
	scripts/utils/customMultiPlot.R\
	scripts/utils/normalizeHic.R\
	scripts/analysis/hoxb13Locus.R
		mkdir -p plots
		Rscript scripts/analysis/hoxb13Locus.R
		
## Same region but zoomed out
plots/hoxb13Locus_v2.pdf:\
	data/diffLoopCounts.rds\
	data/mergedLoopCounts.h5\
	data/raw/hic/genotype/YY1P_22RV1_WT_inter_30.hic\
	data/raw/hic/genotype/YY1P_22RV1_KO_inter_30.hic\
	data/raw/signal/ATAC_seq_EV.bw\
	data/raw/signal/ATAC_seq_YY1KD.bw\
	data/raw/signal/Chip_seq_H3K27ac_.bw\
	data/raw/signal/Chip_seq_INPUT.bw\
	data/raw/signal/Chip_seq_YY1.bw\
	data/raw/signal/RNA_seq_EV.bw\
	data/raw/signal/RNA_seq_YY1KD.bw\
	scripts/utils/customMultiPlot.R\
	scripts/utils/normalizeHic.R\
	scripts/analysis/hoxb13Locus_v2.R
		mkdir -p plots
		Rscript scripts/analysis/hoxb13Locus_v2.R

##############################################
## YY1 enrichment at diff loop anchors      ##
##############################################

## YY1 enrichment at anchors
plots/yy1Enrichment.pdf:\
	tables/diffLoopThresh.txt\
	data/diffLoopCounts.rds\
	data/raw/signal/Chip_seq_YY1.bw\
	scripts/analysis/yy1Enrichment.R
		mkdir -p plots
		Rscript scripts/analysis/yy1Enrichment.R

## CTCF enrichment at anchors (as control)
plots/ctcfEnrichment.pdf:\
	tables/diffLoopThresh.txt\
	data/diffLoopCounts.rds\
	data/raw/signal/CTCF_Chipseq_signal_p-value_GRCh38.bigWig\
	scripts/analysis/ctcfEnrichment.R
		mkdir -p plots
		Rscript scripts/analysis/ctcfEnrichment.R
		
## YY1 enrichment internally (within loop anchors)
plots/yy1InternalEnrichment.pdf:\
	tables/diffLoopThresh.txt\
	data/diffLoopCounts.rds\
	data/raw/signal/Chip_seq_YY1.bw\
	scripts/analysis/yy1InternalEnrichment.R
		mkdir -p plots
		Rscript scripts/analysis/yy1InternalEnrichment.R
		
		
############################################
## How are loops nested? Are gained loops ##
## forming inside lost loops more than    ##
## expected or vice-versa?                ##
############################################

## Permutation plots of loop nesting
plots/loopNesting.pdf:\
	tables/diffLoopThresh.txt\
	data/diffLoopCounts.rds\
	scripts/analysis/loopNesting.R
		mkdir -p plots
		Rscript scripts/analysis/loopNesting.R


###################################
## Deeptools code (not included) ##
###################################

## Export files for deeptools
data/diffLoopBedRegions/static_p0.01_lfc0.bed\
data/diffLoopBedRegions/gained_p0.01_lfc0.bed\
data/diffLoopBedRegions/lost_p0.01_lfc0.bed:\
	data/diffLoopCounts.rds\
	scripts/processing/exportFilesForDeeptools.R
		mkdir -p data/diffLoopBedRegions
		Rscript scripts/processing/exportFilesForDeeptools.R
		
## Create data and heatmap with deeptools
data/deeptools_matrix.mat.gz\
plots/ExampleHeatmap1.png:\
	data/raw/signal/Chip_seq_YY1.bw\
	data/raw/signal/CTCF_Chipseq_signal_p-value_GRCh38.bigWig\
	data/diffLoopBedRegions/static_p0.01_lfc0.bed\
	data/diffLoopBedRegions/gained_p0.01_lfc0.bed\
	data/diffLoopBedRegions/lost_p0.01_lfc0.bed\
	scripts/analysis/deeptools.sh
		mkdir -p data plots
		sh scripts/analysis/deeptools.sh