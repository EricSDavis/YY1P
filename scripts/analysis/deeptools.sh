module load deeptools/3.5.1

# run compute matrix to collect the data needed for plotting
computeMatrix scale-regions \
      -S \
      data/raw/signal/Chip_seq_YY1.bw \
      data/raw/signal/CTCF_Chipseq_signal_p-value_GRCh38.bigWig \
      -R \
      data/diffLoopBedRegions/*.bed \
      --beforeRegionStartLength 5000 \
      --regionBodyLength 10000 \
      --afterRegionStartLength 5000 \
      --skipZeros \
      -o data/deeptools_matrix.mat.gz

plotHeatmap -m data/deeptools_matrix.mat.gz \
      -out plots/ExampleHeatmap1.png \
      