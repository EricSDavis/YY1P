.PHONY: clean

objects :=\
	data/mergedLoops.rds\
	data/mergedLoopCounts.h5\
	data/mergedLoopCounts.rds\
	data/diffLoopCounts.rds\
	plots/surveyDiffLoops.pdf\
	plots/hoxb13Locus.pdf

all: $(objects)

clean:
	rm -rf $(objects)

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

data/diffLoopCounts.rds:\
	data/mergedLoopCounts.h5\
	data/mergedLoopCounts.rds\
	scripts/analysis/differentialLoops.R
		mkdir -p data
		Rscript scripts/analysis/differentialLoops.R

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