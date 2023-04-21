source("scripts/utils/normalizeHic.R")

test_that("normalizeHic works", {
  
  ## Input path of hicFiles
  hicFiles <- c(
    "WT" = "data/raw/hic/genotype/YY1P_22RV1_WT_inter_30.hic",
    "KO" = "data/raw/hic/genotype/YY1P_22RV1_KO_inter_30.hic"
  )

  normalized <- normalizeHic(
    files=hicFiles,
    chrom="chr17",
    chromstart=48702638,
    chromend=48852638
  )
  
  expect_identical(length(normalized), 2L)
  
  
})
