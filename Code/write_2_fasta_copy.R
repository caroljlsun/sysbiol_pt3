library("useful")
library("seqinr")
library("reshape")
library("tidyverse")


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

# Write to fasta files for w/wo G4s
write.fasta(sequences = as.list(seqs_G4_pos$seq), names=rownames(seqs_G4_pos), file.out=paste0(baseDir, "seqs_G4_pos.fasta") )
write.fasta(sequences = as.list(seqs_G4_neg$seq), names=rownames(seqs_G4_neg), file.out=paste0(baseDir, "seqs_G4_neg.fasta") )

# Write a fasta for both postive and negatives


write.fasta(sequences = as.list(seqs_G4_neg$seq), names=rownames(seqs_G4_neg), file.out=paste0(baseDir, "seqs_G4_all.fasta") )