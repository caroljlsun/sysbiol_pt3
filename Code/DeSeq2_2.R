library(DESeq2)

baseDir <- "/home/cjls4/feature_vectors/"
setwd(baseDir)

load("/home/cjls4/feature_vectors/btmf_df.RData")

baseDir <- "/storage/Russell/G4_Tests_Russell"
setwd(baseDir)


DESeqDataSetFromFeatureCounts <- function (sampleTable, directory = ".", design, ignoreRank = FALSE, ...) 
{
  # From https://www.biostars.org/p/277316/
  if (missing(design)) 
    stop("design is missing")
  l <- lapply(as.character(sampleTable[, 2]), function(fn) read.table(file.path(directory, fn), skip=2))
  if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1)))) 
    stop("Gene IDs (first column) differ between files.")
  tbl <- sapply(l, function(a) a$V8) # changes 7 > to 8
  colnames(tbl) <- sampleTable[, 1]
  rownames(tbl) <- l[[1]]$V1
  
  rownames(sampleTable) <- sampleTable[, 1]
  dds <- DESeqDataSetFromMatrix(countData = tbl, colData = sampleTable[, -(1:2), drop = FALSE], design = design, ignoreRank, ...)
  return(dds)
}

#
# Find the featureCount files and then use the filenames to extract the metadata - cell and sex
#
sampleFiles    <- list.files("/storage/CarolSun/Expression", 
                             pattern='*featureCounts.txt$', recursive=T, full.names=T)

sampleNames    <- gsub(".storage.CarolSun.Expression.", "", sampleFiles)
sampleNames    <- gsub(".featureCounts.txt", "", sampleNames)

sampleCellType <- sampleNames
sampleCellType <- gsub("_.*", "", sampleCellType)

sampleSex      <- sampleNames
sampleSex      <- gsub("_R.*", "", sampleSex)
sampleSex      <- gsub(".cell_", "", sampleSex)

sampleTable    <- data.frame(sampleNames=sampleNames, fileNameDGE=sampleFiles, 
                             CellType=sampleCellType, Sex=sampleSex)


#
# Read into DESeq2 and run DE analysis
#


dds        <- DESeqDataSetFromFeatureCounts(sampleTable=sampleTable, directory="", design=~CellType+Sex  )
dds        <- DESeq(dds, parallel=F)

# Use rlog transformed data for a PCA - there are packages to make much better looking PCAs
rld        <- rlog(dds)

plotPCA(rld, intgroup = "CellType")
plotPCA(rld, intgroup = c("CellType", "Sex") )


# compare groups

# Find the groups to be compared based on the input model design
resultsNames(dds)

# basic examples here:
res.celltype <- results(dds, name="CellType_tcell_vs_bcell")
res.sex      <- results(dds, name="Sex_M_vs_F")


plotMA(res.celltype, ylim=c(-5,5))
plotMA(res.sex,     ylim=c(-5,5))


resSig.celltype <- subset(res.celltype, padj < 0.01 & abs(log2FoldChange) >= 2)
head(resSig.celltype)

resSig.sex <- subset(res.sex, padj < 0.01 & abs(log2FoldChange) >= 2)
head(resSig.sex)

#
# I didn't touch your code below this line
#

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
