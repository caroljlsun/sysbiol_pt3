library("useful")
library("seqinr")
library("reshape")
library("tidyverse")
library("M3C") # has the umap functions
library("ggplot2")

baseDir <- "/home/cjls4/kmer_analysis/" # Add where your data/files are


#
# Read in the promoter sequences and write to fasta - also get the annotations for later
#


load("/home/cjls4/feature_vectors/ALL_STRAND_ALL_G4_SEQ.RData")
str(ALL_STRAND_ALL_G4_SEQ)

ALL_STRAND_ALL_G4_SEQ$strand <- gsub("\\-", -1, ALL_STRAND_ALL_G4_SEQ$strand)
ALL_STRAND_ALL_G4_SEQ$strand <- gsub("\\+",  1, ALL_STRAND_ALL_G4_SEQ$strand)
head(ALL_STRAND_ALL_G4_SEQ)

# Extract the positives, write a fasta and annotation df
seqs_G4_pos           <- subset(ALL_STRAND_ALL_G4_SEQ, G4 == 1)
rownames(seqs_G4_pos) <- paste0(seqs_G4_pos$row_names, ".", seqs_G4_pos$G4, ".", seqs_G4_pos$strand)
seqs_G4_pos.annot     <- seqs_G4_pos[, c("G4", "strand")]
head(seqs_G4_pos.annot)

# Extract the negatives, write a fasta and annotation df
seqs_G4_neg           <- subset(ALL_STRAND_ALL_G4_SEQ, G4 == 0)
rownames(seqs_G4_neg) <- paste0(seqs_G4_neg$row_names, ".", seqs_G4_neg$G4, ".", seqs_G4_neg$strand)
seqs_G4_neg.annot     <- seqs_G4_neg[, c("G4", "strand")]
head(seqs_G4_neg.annot)

# Check for duplicate names
intersect( rownames(seqs_G4_pos), rownames(seqs_G4_neg)  )

# Get sequence names linked to their G4 status and strand
seqs_G4_all.annot <- rbind(seqs_G4_pos.annot, seqs_G4_neg.annot)
head(seqs_G4_all.annot)

#
# Need to have run kmer-count and the parser perl script
#
kmer6     <- read.table(paste0(baseDir, "seqs_G4_all.kmer.6.tab"), sep="\t", header=T, row.names = 1)
kmer6.pca <- prcomp(kmer6, scale. = TRUE, center = TRUE)
corner(kmer6.pca$x)

pca.data.G4           <- merge(as.data.frame(kmer6.pca$x)[, c(1:10)], seqs_G4_all.annot, by=0)
rownames(pca.data.G4) <- pca.data.G4$Row.names
pca.data.G4           <- pca.data.G4[, c(2:ncol(pca.data.G4))]

ggplot(pca.data.G4, aes(x=PC1, y=PC2) ) +
  geom_point(aes(colour=as.factor(G4)), size=2.5, alpha=.1, shape=20) +
  scale_colour_viridis_d(alpha = 0.1, option = "D",
                         begin = 0.1, end = 0.7, 
                         name = "G4", labels = c("negative", "positive"))+
  #scale_colour_manual( name="G4", values=c("0"="darkgreen", "1"="magenta") ) +
  ggtitle("PCA of k-mers") +
  theme_bw()



#
# Originally formatted like this for another program, but work for the UMAP input
#
kmer6Data <- list(data=kmer6, G4=seqs_G4_all.annot$G4, strand=seqs_G4_all.annot$strand, label=paste0(seqs_G4_all.annot$G4, "_", seqs_G4_all.annot$strand) )
kmer6Data$label <- gsub("1_1","G4pos_plus",    kmer6Data$label)
kmer6Data$label <- gsub("1_-1","G4pos_minus",  kmer6Data$label)
kmer6Data$label <- gsub("0_1","G4pneg_plus",   kmer6Data$label)
kmer6Data$label <- gsub("0_-1","G4neg_minus",  kmer6Data$label)
corner(kmer6Data$data)

#kmer6Data$

umap   <- umap::umap(as.matrix(kmer6Data$data))
scores <- data.frame(umap$layout)  
head(scores)


# hospital_names <- c(
#   `Hospital#1` = "Some Hospital",
#   `Hospital#2` = "Another Hospital",
#   `Hospital#3` = "Hospital Number 3",
#   `Hospital#4` = "The Other Hospital"
# )

# facet_grid(hospital ~ ., labeller = as_labeller(hospital_names))

umap_names <- c(`0` = "G4 negative",
                `1` = "G4 positive")



uh <- ggplot(data = scores, aes(x = X1, y = X2, label = kmer6Data$G4) ) + 
  geom_point(aes(colour = as.factor(kmer6Data$label)), size = 1, alpha=.5) + 
  scale_colour_viridis_d(alpha = 0.07, option = "D",
                         direction = -1, begin = 0,
                         end = 0.4, name = "Type of K-mer", 
                         labels = c("G4 negative, minus strand",
                                    "G4 negative, plus strand",
                                    "G4 positive, minus strand", 
                                    "G4 positive, plus strand"))+
  facet_wrap( ~ kmer6Data$G4, labeller = as_labeller(umap_names) ) +
  theme_bw()+
  ggtitle("K-mer analysis of promoter similarity")+
  xlab("UMAP1")+
  ylab("UMAP2")

uh

save(umap, file = "umap.RData")

