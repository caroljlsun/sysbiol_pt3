##### G quadruplex data #####

load("/storage/CarolSun/G4s/plus_PDS_workspace_g4.RData")

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

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

library(DNAshapeR)
library(seqinr)
library(fields)
library(Biostrings)


##### Global options #####

setwd(dir = "/storage/CarolSun/G4s/")
options(scipen=999)
options(expressions = 500000)

##### Importing data #####

mouse_plus_PDS <- readGeneric("Mouse_all_w15_th-1_plus.hits.max.PDS.w50.35.bed",
                            chr = 1,
                            start = 2,
                            end = 3,
                            strand = NULL)

strand(mouse_plus_PDS) <- "-"


# experimentally observed Qs reverse strand in PDS

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

intersect_plus_PDS <- findOverlapPairs(mouse_plus_PDS, mouse_genes_promoters_coords, type = "within")

#check dimensions

head(intersect_plus_PDS)
head(intersect_plus_PDS@first)
head(intersect_plus_PDS@second)
length(intersect_plus_PDS@first) # 51022

length(mouse_genes_promoters_coords)#63562
length(mouse_plus_PDS)#873853

#get sequences

sequences_plus_PDS <- getSeq(mouse_genome, intersect_plus_PDS@second)
head(sequences_plus_PDS)

#get promoter coordinates
head(intersect_plus_PDS@second@ranges)

promoter_plus_PDS_start <- intersect_plus_PDS@second@ranges@start
promoter_plus_PDS_end <- intersect_plus_PDS@second@ranges@start + intersect_plus_PDS@second@ranges@width - 1
# check width

plus_PDS_promoter_widths <- promoter_plus_PDS_end - promoter_plus_PDS_start + 1# surely +1 right

#get G4 coordinates

g_plus_PDS_start <- intersect_plus_PDS@first@ranges@start
g_plus_PDS_end <- intersect_plus_PDS@first@ranges@start + intersect_plus_PDS@first@ranges@width - 1


# make contextual G4 coordinates
# so that we get the start and ends of the G4s in the context of the promoters themselves
# with the beginning being the start of the promoter etc.

context_g_plus_PDS_start <- g_plus_PDS_start - promoter_plus_PDS_start
context_g_plus_PDS_end <- g_plus_PDS_end - promoter_plus_PDS_start

#make a featurevector

# first get some IDs for the promoters

plus_PDS_promoter_chromosome <- intersect_plus_PDS@second@elementMetadata$tx_chrom

plus_PDS_promoter_coord_names <- paste0(plus_PDS_promoter_chromosome, "_", promoter_plus_PDS_start, "_", promoter_plus_PDS_end, collapse = ";")

plus_PDS_promoter_names <- as.list(unlist(strsplit(plus_PDS_promoter_coord_names, "[;]")))

# they're not unique though, so may need to use tx_ids as well

plus_PDS_unique_names <- intersect_plus_PDS@second@elementMetadata$tx_name

##### assignment of specific G4 positions to each promoter #####

plus_PDS_sequences_list <- as.character(plus_PDS_promoter_names, use.names=TRUE)
plus_PDS_sequence_names <- as_tibble(plus_PDS_sequences_list)

plus_PDS_tibble <- plus_PDS_sequence_names %>% 
  add_column(context_g_plus_PDS_start) %>% 
  add_column(context_g_plus_PDS_end) %>% 
  add_column(as.character(sequences_plus_PDS)) %>% 
  add_column(plus_PDS_unique_names)%>% 
  dplyr::rename(chromosome = value) %>% 
  dplyr::rename(g_start = context_g_plus_PDS_start) %>% 
  dplyr::rename(g_end = context_g_plus_PDS_end) %>% 
  dplyr::rename(seq = `as.character(sequences_plus_PDS)`)

##### featurevector #####

plus_PDS_fv <- matrix( data = 0, nrow = length(context_g_plus_PDS_start), 
                     ncol = 2750)

for (i in 1:dim(plus_PDS_fv)[1]){
  plus_PDS_fv[i, (plus_PDS_tibble$g_start[i]:plus_PDS_tibble$g_end[i])] <- 1
}

#check

table(plus_PDS_fv[1,])
table(plus_PDS_fv[90,])

#convert to data.table

plus_PDS_fv <- as.data.table(plus_PDS_fv)


# get the tibble, since you need to match the sequences etc.

plus_PDS_fv_tibble <- plus_PDS_fv %>% 
  add_column(plus_PDS_tibble$chromosome, .before = 1) %>% 
  add_column(plus_PDS_tibble$seq, .before = 2) %>% 
  dplyr::rename(seq = `plus_PDS_tibble$seq`) %>% 
  dplyr::rename(chromosome = `plus_PDS_tibble$chromosome`)


#check for duplicates...
length(unique(plus_PDS_fv_tibble$seq)) #16161

#right, so there are both duplicates AND multiple instances of G4s that haven't merged together
#will need to use the coalesce stuff again lol

#for all chromosomes?

plus_PDS_seq <- tibble(plus_PDS_fv_tibble$seq) %>% 
  arrange(plus_PDS_fv_tibble$seq) %>% 
  group_by(plus_PDS_fv_tibble$seq)

length(unique(plus_PDS_fv_tibble$seq))

plus_PDS_numbers <- tibble(plus_PDS_fv_tibble) %>% 
  arrange(seq)

plus_PDS_numbers <- plus_PDS_numbers %>% 
  dplyr::select(-chromosome, -seq)

indices_pds <- group_rows(plus_PDS_seq)
coalesced_pds <- vector(length = 1)


plus_PDS_numbers[plus_PDS_numbers == 0] <- NA
plus_PDS_numbers <- lapply(plus_PDS_numbers, as.numeric)
plus_PDS_numbers <- transpose(plus_PDS_numbers)

head(indices_pds)
does_it_work <- fcoalesce(plus_PDS_numbers[c(10334, 10335, 10336, 10337, 10338, 10339)])

#does_it_work <- fcoalesce(plus_PDS_numbers[10334], plus_PDS_numbers[10335])

table(does_it_work)
table(plus_PDS_numbers[10334])

istart <- lapply(indices_pds, min)
ihold <- lapply(indices_pds, length)
iend <- mapply(function(x,y)list(sum(x, y, -1)),x = istart, y = ihold)

istart[[1]]:iend[[1]]

pds_coalesce <- function(i) {
  coalesced_pds <- list(as.numeric(c(fcoalesce(plus_PDS_numbers[c(indices_pds[[i]])]))))
  coalesced_pds
}

pds <- sapply(1:length(indices_pds), pds_coalesce)

pds_matrix <- ldply(pds)

pds_matrix[is.na(pds_matrix)] <- 0


# will need to somehow make one that is a simple absence presence vector
# for the use of other feature vectors

sorted_chr_seq <- plus_PDS_fv_tibble %>% 
  arrange(plus_PDS_fv_tibble$seq) %>% 
  distinct(seq, .keep_all = T)

sorted_chr_seq <- sorted_chr_seq %>% 
  dplyr::select(chromosome, seq)

plus_PDS_fv_final <- pds_matrix %>% 
  add_column(sorted_chr_seq, .before = 1)

save.image("plus_PDS_workspace_g4.RData")

setwd("/home/cjls4/feature_vectors/")

save(plus_PDS_fv_final, file = "plus_PDS_G4_final")

