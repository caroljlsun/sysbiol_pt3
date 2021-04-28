##### G quadruplex data #####

load("/storage/CarolSun/G4s/minus_PDS_workspace_g4.RData")

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

mouse_minus_PDS <- readGeneric("Mouse_all_w15_th-1_minus.hits.max.PDS.w50.35.bed",
                             chr = 1,
                             start = 2,
                             end = 3,
                             strand = NULL)

strand(mouse_minus_PDS) <- "+"


# experimentally observed Qs forward strand in PDS

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

intersect_minus_PDS <- findOverlapPairs(mouse_minus_PDS, mouse_genes_promoters_coords, type = "within")

#check dimensions

head(intersect_minus_PDS)
head(intersect_minus_PDS@first)
head(intersect_minus_PDS@second)
length(intersect_minus_PDS@first) # 51956

length(mouse_genes_promoters_coords)#63562
length(mouse_minus_PDS)#873010

#get sequences

sequences_minus_PDS <- getSeq(mouse_genome, intersect_minus_PDS@second)
head(sequences_minus_PDS)

#get promoter coordinates
head(intersect_minus_PDS@second@ranges)

promoter_minus_PDS_start <- intersect_minus_PDS@second@ranges@start
promoter_minus_PDS_end <- intersect_minus_PDS@second@ranges@start + intersect_minus_PDS@second@ranges@width - 1
# check width

minus_PDS_promoter_widths <- promoter_minus_PDS_end - promoter_minus_PDS_start + 1# surely +1 right

#get G4 coordinates

g_minus_PDS_start <- intersect_minus_PDS@first@ranges@start
g_minus_PDS_end <- intersect_minus_PDS@first@ranges@start + intersect_minus_PDS@first@ranges@width - 1


# make contextual G4 coordinates
# so that we get the start and ends of the G4s in the context of the promoters themselves
# with the beginning being the start of the promoter etc.

context_g_minus_PDS_start <- g_minus_PDS_start - promoter_minus_PDS_start
context_g_minus_PDS_end <- g_minus_PDS_end - promoter_minus_PDS_start

#make a featurevector

# first get some IDs for the promoters

minus_promoter_chromosome <- intersect_minus_PDS@second@elementMetadata$tx_chrom

minus_PDS_promoter_coord_names <- paste0(minus_promoter_chromosome, "_", promoter_minus_PDS_start, "_", promoter_minus_PDS_end, collapse = ";")

minus_PDS_promoter_names <- as.list(unlist(strsplit(minus_PDS_promoter_coord_names, "[;]")))

# they're not unique though, so may need to use tx_ids as well

minus_PDS_unique_names <- intersect_minus_PDS@second@elementMetadata$tx_name

##### assignment of specific G4 positions to each promoter #####

minus_PDS_sequences_list <- as.character(minus_PDS_promoter_names, use.names=TRUE)
minus_PDS_sequence_names <- as_tibble(minus_PDS_sequences_list)

minus_PDS_tibble <- minus_PDS_sequence_names %>% 
  add_column(context_g_minus_PDS_start) %>% 
  add_column(context_g_minus_PDS_end) %>% 
  add_column(as.character(sequences_minus_PDS)) %>% 
  add_column(minus_PDS_unique_names)%>% 
  dplyr::rename(chromosome = value) %>% 
  dplyr::rename(g_start = context_g_minus_PDS_start) %>% 
  dplyr::rename(g_end = context_g_minus_PDS_end) %>% 
  dplyr::rename(seq = `as.character(sequences_minus_PDS)`)

##### featurevector #####

minus_PDS_fv <- matrix( data = 0, nrow = length(context_g_minus_PDS_start), 
                      ncol = 2750)

for (i in 1:dim(minus_PDS_fv)[1]){
  minus_PDS_fv[i, (minus_PDS_tibble$g_start[i]:minus_PDS_tibble$g_end[i])] <- 1
}

#check

table(minus_PDS_fv[1,])
table(minus_PDS_fv[90,])

#convert to data.table

minus_PDS_fv <- as.data.table(minus_PDS_fv)


# get the tibble, since you need to match the sequences etc.

minus_PDS_fv_tibble <- minus_PDS_fv %>% 
  add_column(minus_PDS_tibble$chromosome, .before = 1) %>% 
  add_column(minus_PDS_tibble$seq, .before = 2) %>% 
  dplyr::rename(seq = `minus_PDS_tibble$seq`) %>% 
  dplyr::rename(chromosome = `minus_PDS_tibble$chromosome`)


#check for duplicates...
length(unique(minus_PDS_fv_tibble$seq)) #16487

#right, so there are both duplicates AND multiple instances of G4s that haven't merged together
#will need to use the coalesce stuff again lol

#for all chromosomes?

minus_PDS_seq <- tibble(minus_PDS_fv_tibble$seq) %>% 
  arrange(minus_PDS_fv_tibble$seq) %>% 
  group_by(minus_PDS_fv_tibble$seq)

length(unique(minus_PDS_fv_tibble$seq))

minus_PDS_numbers <- tibble(minus_PDS_fv_tibble) %>% 
  arrange(seq)

minus_PDS_numbers <- minus_PDS_numbers %>% 
  dplyr::select(-chromosome, -seq)

indices_pds <- group_rows(minus_PDS_seq)
coalesced_pds <- vector(length = 1)


minus_PDS_numbers[minus_PDS_numbers == 0] <- NA
minus_PDS_numbers <- lapply(minus_PDS_numbers, as.numeric)
minus_PDS_numbers <- transpose(minus_PDS_numbers)

head(indices_pds)
does_it_work <- fcoalesce(minus_PDS_numbers[c(6411,6422,6428)])

istart <- lapply(indices_pds, min)
ihold <- lapply(indices_pds, length)
iend <- mapply(function(x,y)list(sum(x, y, -1)),x = istart, y = ihold)

istart[[1]]:iend[[1]]

pds_coalesce <- function(i) {
  coalesced_pds <- list(as.numeric(c(fcoalesce(minus_PDS_numbers[c(indices_pds[[i]])]))))
  coalesced_pds
}

pds <- sapply(1:length(indices_pds), pds_coalesce)

pds_matrix <- ldply(pds)

pds_matrix[is.na(pds_matrix)] <- 0


# will need to somehow make one that is a simple absence presence vector
# for the use of other feature vectors

sorted_chr_seq <- minus_PDS_fv_tibble %>% 
  arrange(minus_PDS_fv_tibble$seq) %>% 
  distinct(seq, .keep_all = T)

sorted_chr_seq <- sorted_chr_seq %>% 
  dplyr::select(chromosome, seq)

minus_PDS_fv_final <- pds_matrix %>% 
  add_column(sorted_chr_seq, .before = 1)

save.image("minus_PDS_workspace_g4.RData")

setwd("/home/cjls4/feature_vectors/")

save(minus_PDS_fv_final, file = "minus_PDS_G4_final")

save(minus_PDS_fv, file = "minus_PDS_G4_dataframe.Rdata")
