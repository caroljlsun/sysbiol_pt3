library(methylKit)
library(GenomicRanges)
library(genomation)
library(readr)
library(BiSeq)

setwd(dir = "/storage/CarolSun/Methylation/")


fem_B_10_10 <- read_tsv("B6_F_B_BS.10.methy.combined_strand.mincov_10.tsv")
fem_B_10_10_conv <- convert("B6_F_B_BS.10.methy.combined_strand.mincov_10.tsv",
                            "B6_F_B_BS.10.methy.combined_strand.mincov_10.txt")

as(fem_B_10_10, "methylRaw")


file <- system.file("extdata", "CpG_context_test_sample.cov", package = "BiSeq")
rrbs <- readBismark(file,
                    colData= DataFrame(row.names="sample_1"))

please <- readBismark("B6_F_B_BS.10.methy.combined_strand.mincov_10.tsv",
            colData= DataFrame(row.names="sample_1"))

#testfiles <- list( "B6_F_B")

tsv_columns <- list()

testDB <- methRead("B6_F_B_BS.10.methy.combined_strand.mincov_10.tsv",
                   sample.id = "FB1010",
                   assembly="black6",
                   dbtype = "tabix",
                   pipeline = "amp"
                   sep = "\t",
                   context ="CpG",
                   dbdir = "methylDB")
methRead(
  location,
  sample.id,
  assembly,
  dbtype,
  pipeline,
  header,
  skip,
  sep,
  context,
  resolution,
  dbdir,
  mincov
)
# A tibble: 6 x 4
# chr1  `3003379`   `0`  `10`
# <chr>     <dbl> <dbl> <dbl>
#   1 chr1    3003582     0    12
# 2 chr1    3003885     0    11
# 3 chr1    3003898     0    11
# 4 chr1    3007532     2     8
# 5 chr1    3007580     1     9
# 6 chr1    3010645     3    11
