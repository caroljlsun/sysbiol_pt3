library(DESeq2)

load("btmf_df.RData")

top_100_names[60] <- "Sfmbt2_1"
top_100_names[90] <- "Sfmbt2_2"

rownames(btmf_df) <- unlist(top_100_names)

coldata <- matrix(data = NA, nrow = 4, ncol = 2)
coldata[,1] <- c("b","b","t","t")
coldata[,2] <- c("f","m","f","m")

# coldata[,3] <- c(0,1,0,1)
# coldata[,4] <- c(1,0,1,0)

rownames(coldata) <- c("bf", "bm", "tf", "tm")
#colnames(coldata) <- c("b_cell", "t_cell", "m", "f")

colnames(coldata) <- c("cell", "sex")

dds_full <- DESeqDataSetFromMatrix(
  countData = btmf_df,
  colData = coldata,
  design = ~ cell + sex)

dds_full <- DESeq(dds_full)

dds_full <- estimateSizeFactors(dds_full)

dds_full <- estimateDispersionsGeneEst(dds_full)

dispersions(dds_full) <- mcols(dds_full)$dispGeneEst

dds_full <- nbinomLRT(dds_full, full = design(dds_full), reduced = ~sex)

res <- results(dds_full)

hmm <- DESeqTransform(dds_full)
plotPCA(hmm, intgroup = "cell")

plotMA(res)
plotDispEsts(dds_full)
hist( res$pvalue, breaks=20, col="grey" )


save(dds_full, file = "dds_full.RData")
save(res, file = "res.RData")
