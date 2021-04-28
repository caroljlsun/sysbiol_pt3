##### G quadruplex data #####

load("/storage/CarolSun/G4s/plus_k_workspace_g4.RData")

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

mouse_plus_k <- readGeneric("Mouse_all_w15_th-1_plus.hits.max.K.w50.25.bed",
                               chr = 1,
                               start = 2,
                               end = 3,
                               strand = NULL)

strand(mouse_plus_k) <- "-"


# experimentally observed Qs reverse strand in K+

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

intersect_plus_k <- findOverlapPairs(mouse_plus_k, mouse_genes_promoters_coords, type = "within")

#check dimensions

head(intersect_plus_k)
head(intersect_plus_k@first)
head(intersect_plus_k@second)
length(intersect_plus_k@first) # 23559

length(mouse_genes_promoters_coords)#63562
length(mouse_plus_k)#398623

#get sequences

sequences_plus_k <- getSeq(mouse_genome, intersect_plus_k@second)
head(sequences_plus_k)

#get promoter coordinates
head(intersect_plus_k@second@ranges)

promoter_plus_k_start <- intersect_plus_k@second@ranges@start
promoter_plus_k_end <- intersect_plus_k@second@ranges@start + intersect_plus_k@second@ranges@width - 1
# check width

plus_k_promoter_widths <- promoter_plus_k_end - promoter_plus_k_start + 1# surely +1 right

#get G4 coordinates

g_plus_k_start <- intersect_plus_k@first@ranges@start
g_plus_k_end <- intersect_plus_k@first@ranges@start + intersect_plus_k@first@ranges@width - 1


# make contextual G4 coordinates
# so that we get the start and ends of the G4s in the context of the promoters themselves
# with the beginning being the start of the promoter etc.

context_g_plus_k_start <- g_plus_k_start - promoter_plus_k_start
context_g_plus_k_end <- g_plus_k_end - promoter_plus_k_start

#make a featurevector

# first get some IDs for the promoters

plus_k_promoter_chromosome <- intersect_plus_k@second@elementMetadata$tx_chrom

plus_k_promoter_coord_names <- paste0(plus_k_promoter_chromosome, "_", promoter_plus_k_start, "_", promoter_plus_k_end, collapse = ";")

plus_k_promoter_names <- as.list(unlist(strsplit(plus_k_promoter_coord_names, "[;]")))

# they're not unique though, so may need to use tx_ids as well

plus_k_unique_names <- intersect_plus_k@second@elementMetadata$tx_name

##### assignment of specific G4 positions to each promoter #####

plus_k_sequences_list <- as.character(plus_k_promoter_names, use.names=TRUE)
plus_k_sequence_names <- as_tibble(plus_k_sequences_list)

plus_k_tibble <- plus_k_sequence_names %>% 
  add_column(context_g_plus_k_start) %>% 
  add_column(context_g_plus_k_end) %>% 
  add_column(as.character(sequences_plus_k)) %>% 
  add_column(plus_k_unique_names)%>% 
  dplyr::rename(chromosome = value) %>% 
  dplyr::rename(g_start = context_g_plus_k_start) %>% 
  dplyr::rename(g_end = context_g_plus_k_end) %>% 
  dplyr::rename(seq = `as.character(sequences_plus_k)`)

##### featurevector #####

plus_k_fv <- matrix( data = 0, nrow = length(context_g_plus_k_start), 
                        ncol = 2750)

for (i in 1:dim(plus_k_fv)[1]){
  plus_k_fv[i, (plus_k_tibble$g_start[i]:plus_k_tibble$g_end[i])] <- 1
}

#check

table(plus_k_fv[1,])
table(plus_k_fv[90,])

#convert to data.table

plus_k_fv <- as.data.table(plus_k_fv)


# get the tibble, since you need to match the sequences etc.

plus_k_fv_tibble <- plus_k_fv %>% 
  add_column(plus_k_tibble$chromosome, .before = 1) %>% 
  add_column(plus_k_tibble$seq, .before = 2) %>% 
  dplyr::rename(seq = `plus_k_tibble$seq`) %>% 
  dplyr::rename(chromosome = `plus_k_tibble$chromosome`)


#check for duplicates...
length(unique(plus_k_fv_tibble$seq)) #10336

#right, so there are both duplicates AND multiple instances of G4s that haven't merged together
#will need to use the coalesce stuff again lol

#for all chromosomes?

plus_k_seq <- tibble(plus_k_fv_tibble$seq) %>% 
  arrange(plus_k_fv_tibble$seq) %>% 
  group_by(plus_k_fv_tibble$seq)

length(unique(plus_k_fv_tibble$seq))

plus_k_numbers <- tibble(plus_k_fv_tibble) %>% 
  arrange(seq)

plus_k_numbers <- plus_k_numbers %>% 
  dplyr::select(-chromosome, -seq)

indices_pk <- group_rows(plus_k_seq)
coalesced_pk <- vector(length = 1)


plus_k_numbers[plus_k_numbers == 0] <- NA
plus_k_numbers <- lapply(plus_k_numbers, as.numeric)
plus_k_numbers <- transpose(plus_k_numbers)

head(indices_pk)
does_it_work <- fcoalesce(plus_k_numbers[c(4854,4855)])
table(does_it_work)
table(plus_k_numbers[4854])

istart <- lapply(indices_pk, min)
ihold <- lapply(indices_pk, length)
iend <- mapply(function(x,y)list(sum(x, y, -1)),x = istart, y = ihold)

istart[[1]]:iend[[1]]

pk_coalesce <- function(i) {
  coalesced_pk <- list(as.numeric(c(fcoalesce(plus_k_numbers[c(indices_pk[[i]])]))))
  coalesced_pk
}

pk <- sapply(1:length(indices_pk), pk_coalesce)

pk_matrix <- ldply(pk)

pk_matrix[is.na(pk_matrix)] <- 0


# will need to somehow make one that is a simple absence presence vector
# for the use of other feature vectors

sorted_chr_seq <- plus_k_fv_tibble %>% 
  arrange(plus_k_fv_tibble$seq) %>% 
  distinct(seq, .keep_all = T)

sorted_chr_seq <- sorted_chr_seq %>% 
  dplyr::select(chromosome, seq)

plus_k_fv_final <- pk_matrix %>% 
  add_column(sorted_chr_seq, .before = 1)

save.image("plus_k_workspace_g4.RData")

setwd("/home/cjls4/feature_vectors/")

save(plus_k_fv_final, file = "plus_k_G4_final")

