.PHONY: clean

objects :=\
	data/mergedLoops.rds

all: $(objects)

clean:
	rm -rf $(objects)

data/mergedLoops.rds:\
	data/raw/loops/KO_inter_30_5KbLoops.txt\
	data/raw/loops/WT_inter_30_5KbLoops.txt\
	data/raw/loops/YY1P_22RV1_KO_1_1_inter_30_5KbLoops.txt\
	data/raw/loops/YY1P_22RV1_KO_1_2_inter_30_5KbLoops.txt\
	data/raw/loops/YY1P_22RV1_KO_2_1_inter_30_5KbLoops.txt\
	data/raw/loops/YY1P_22RV1_KO_2_2_inter_30_5KbLoops.txt\
	data/raw/loops/YY1P_22RV1_WT_1_1_inter_30_5KbLoops.txt\
	data/raw/loops/YY1P_22RV1_WT_1_2_inter_30_5KbLoops.txt\
	data/raw/loops/YY1P_22RV1_WT_2_1_inter_30_5KbLoops.txt\
	data/raw/loops/YY1P_22RV1_WT_2_2_inter_30_5KbLoops.txt\
	scripts/processing/mergeLoops.R
		mkdir -p data
		Rscript scripts/processing/mergeLoops.R 