##### G quadruplex data #####

#general idea is to get a 4 choice one hot encode
#no G, G with K, G with Li, G with both K and Li

load("/storage/CarolSun/G4s/minus_k_workspace_g4.RData")

##### Libraries #####
library(GenomicFeatures)
library(genomation)
library(GenomicRanges)
library(tidyverse)
library(tibble)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(useful)
library(data.table)
library(plyr)
##### Global options #####

setwd(dir = "/storage/CarolSun/G4s/")
options(scipen=999)
options(expressions = 500000)

##### Importing data #####

mouse_minus_K <- readGeneric("Mouse_all_w15_th-1_minus.hits.max.K.w50.25.bed",
                             chr = 1,
                             start = 2,
                             end = 3,
                             strand = NULL)

strand(mouse_minus_K) <- "+"


# experimentally observed Qs forward strand in K
# maximal percentage mis- match value above 25 for K+ and above 35 for K++PDS (all the max files)

# check their min coverage requirements!

##### Intersect #####

mouse_annotation <- TxDb.Mmusculus.UCSC.mm10.knownGene
mouse_genome <- BSgenome.Mmusculus.UCSC.mm10

seqlevels(mouse_annotation) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                 "chr11","chr12","chr13","chr14","chr15", "chr16","chr17","chr18","chr19","chrX","chrY")

mouse_genes <- keys(mouse_annotation)

mouse_genes_transcript_coords <- transcriptsBy(mouse_annotation, by = "gene")

mouse_genes_promoters <- getPromoterSeq(mouse_genes_transcript_coords, mouse_genome, upstream = 2500, downstream = 250)

mouse_genes_promoters_coords   <- promoters(mouse_annotation, upstream=2500, downstream=250, 
                                            columns=c("tx_id", "tx_chrom", "tx_strand", "tx_name", "gene_id"))

intersect_minus_k <- findOverlapPairs(mouse_minus_K, mouse_genes_promoters_coords, type = "within")

#check dimensions

head(intersect_minus_k)
head(intersect_minus_k@first)
head(intersect_minus_k@second)
length(intersect_minus_k@first) # 23718
length(intersect_minus_k@second)#23718

length(mouse_genes_promoters_coords)#63562
length(mouse_minus_K)#399166

#get sequences

sequences_minus_k <- getSeq(mouse_genome, intersect_minus_k@second)
head(sequences_minus_k)

#get promoter coordinates
head(intersect_minus_k@second@ranges)

promoter_minus_k_start <- intersect_minus_k@second@ranges@start
promoter_minus_k_end <- intersect_minus_k@second@ranges@start + intersect_minus_k@second@ranges@width - 1
# check width

minus_k_promoter_widths <- promoter_minus_k_end - promoter_minus_k_start + 1# surely +1 right

#get G4 coordinates

g_minus_k_start <- intersect_minus_k@first@ranges@start
g_minus_k_end <- intersect_minus_k@first@ranges@start + intersect_minus_k@first@ranges@width - 1


# make contextual G4 coordinates
# so that we get the start and ends of the G4s in the context of the promoters themselves
# with the beginning being the start of the promoter etc.

context_g_minus_k_start <- g_minus_k_start - promoter_minus_k_start
context_g_minus_k_end <- g_minus_k_end - promoter_minus_k_start

#make a featurevector

# first get some IDs for the promoters

minus_promoter_chromosome <- intersect_minus_k@second@elementMetadata$tx_chrom

minus_k_promoter_coord_names <- paste0(minus_promoter_chromosome, "_", promoter_minus_k_start, "_", promoter_minus_k_end, collapse = ";")

minus_k_promoter_names <- as.list(unlist(strsplit(minus_k_promoter_coord_names, "[;]")))

# they're not unique though, so may need to use tx_ids as well

minus_k_unique_names <- intersect_minus_k@second@elementMetadata$tx_name

##### assignment of specific G4 positions to each promoter #####

minus_k_sequences_list <- as.character(minus_k_promoter_names, use.names=TRUE)
minus_k_sequence_names <- as_tibble(minus_k_sequences_list)

minus_k_tibble <- minus_k_sequence_names %>% 
  add_column(context_g_minus_k_start) %>% 
  add_column(context_g_minus_k_end) %>% 
  add_column(as.character(sequences_minus_k)) %>% 
  add_column(minus_k_unique_names) %>% 
  rename(chromosome = value) %>% 
  rename(g_start = context_g_minus_k_start) %>% 
  rename(g_end = context_g_minus_k_end) %>% 
  rename(seq = `as.character(sequences_minus_k)`)

##### featurevector #####

minus_k_fv <- matrix( data = 0, nrow = length(context_g_minus_k_start), 
                         ncol = 2750)

for (i in 1:dim(minus_k_fv)[1]){
  minus_k_fv[i, (minus_k_tibble$g_start[i]:minus_k_tibble$g_end[i])] <- 1
}

#check

table(minus_k_fv[1,])
table(minus_k_fv[600,])

#convert to data.table

minus_k_fv <- as.data.table(minus_k_fv)


#HOLD ON, MAKE IT DROP
#lol
# get the tibble, since you need to match the sequences etc.

#check for duplicates...

length(unique(minus_k_fv_tibble$seq))
qc <- minus_k_fv_tibble[1:10,]
qc_1 <- minus_k_tibble[1:10,]
#right, so there are both duplicates AND multiple instances of G4s that haven't merged together
#will need to use the coalesce stuff again lol

please_qc <- tibble(qc$seq) %>% 
  arrange(qc$seq) %>% 
  group_by(qc$seq)

qc_numbers <- tibble(qc) %>% 
  arrange(seq)

qc_numbers <- qc_numbers %>% 
  dplyr::select(-chromosome, -seq)

indices_qc <- group_rows(please_qc)
coalesced_qc <- vector(length = 1)


qc_numbers[qc_numbers == 0] <- NA
qc_numbers <- lapply(qc_numbers, as.numeric)
qc_numbers <- transpose(qc_numbers)



does_it_work <- fcoalesce(qc_numbers[c(2,3)])

istart <- lapply(indices_qc, min)
ihold <- lapply(indices_qc, length)
iend <- mapply(function(x,y)list(sum(x, y, -1)),x = istart, y = ihold)

istart[[1]]:iend[[1]]

qc_coalesce <- function(i) {
  coalesced_qc <- list(as.numeric(c(fcoalesce(qc_numbers[istart[[i]]:iend[[i]]]))))
  coalesced_qc
}

qc_3 <- sapply(1:length(indices_qc), qc_coalesce)

c3_matrix <- ldply(qc_3)

c3_matrix[is.na(c3_matrix)] <- 0

#repeat the above for all chromosomes?

minus_k_seq <- tibble(minus_k_fv_tibble$seq) %>% 
  arrange(minus_k_fv_tibble$seq) %>% 
  group_by(minus_k_fv_tibble$seq)

length(unique(minus_k_fv_tibble$seq))

minus_k_numbers <- tibble(minus_k_fv_tibble) %>% 
  arrange(seq)

minus_k_numbers <- minus_k_numbers %>% 
  dplyr::select(-chromosome, -seq)

indices_mk <- group_rows(minus_k_seq)
coalesced_mk <- vector(length = 1)


minus_k_numbers[minus_k_numbers == 0] <- NA
minus_k_numbers <- lapply(minus_k_numbers, as.numeric)
minus_k_numbers <- transpose(minus_k_numbers)

does_it_work <- fcoalesce(minus_k_numbers[c(2942,2954)])

istart <- lapply(indices_mk, min)
ihold <- lapply(indices_mk, length)
iend <- mapply(function(x,y)list(sum(x, y, -1)),x = istart, y = ihold)

istart[[1]]:iend[[1]]

mk_coalesce <- function(i) {
  coalesced_mk <- list(as.numeric(c(fcoalesce(minus_k_numbers[c(indices_mk[[i]])]))))
  coalesced_mk
}

mk <- sapply(1:length(indices_mk), mk_coalesce)

mk_matrix <- ldply(mk)

mk_matrix[is.na(mk_matrix)] <- 0


# will need to somehow make one that is a simple absence presence vector
# for the use of other feature vectors

sorted_chr_seq <- minus_k_fv_tibble %>% 
  arrange(minus_k_fv_tibble$seq) %>% 
  distinct(seq, .keep_all = T)

sorted_chr_seq <- sorted_chr_seq %>% 
  dplyr::select(chromosome, seq)

#minus_k_fv_tibble <- mk_matrix %>% 
  #add_column(minus_k_tibble$chromosome, .before = 1) %>% 
  #add_column(minus_k_tibble$seq, .before = 2) %>% 
  #dplyr::rename(seq = `minus_k_tibble$seq`) %>% 
  #dplyr::rename(chromosome = `minus_k_tibble$chromosome`)


minus_k_fv_final <- mk_matrix %>% 
  add_column(sorted_chr_seq, .before = 1)

save.image("minus_k_workspace_g4.RData")

setwd("/home/cjls4/feature_vectors/")

save(minus_k_fv_final, file = "minus_k_G4_final")
save(minus_k_fv, file = "minus_k_G4_dataframe.Rdata")
