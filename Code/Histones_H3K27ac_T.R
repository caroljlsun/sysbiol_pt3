# if you want the original workspace when I was working on this please use ls -a to find a file called .RData
# if it exists, and it should, use something like load("/storage/CarolSun/ChIP/H3K27ac_T.RData")

load("H3K27ac_T.RData")

library(GenomicFeatures)
library(genomation)
library(GenomicRanges)
library(tidyverse)
library(tibble)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(useful)
library(plyr)

##### Global options #####

setwd("/storage/CarolSun/ChIP/")
options(scipen=999)
options(expressions = 500000)
#browseURL('https://www.youtube.com/watch?v=QH2-TGUlwu4')

##### Import #####

H3K27ac_T <- read_bigwig("GSE94658_H3K27ac_NaiveT_BL6_dedup_merged.bw")

##### Quick inspection #####
head(H3K27ac_T)
max(H3K27ac_T$score) #13.87472
min(H3K27ac_T$score) #0

hist(H3K27ac_T$score)

table(H3K27ac_T@strand) #absolutely no strand information, which is to be expected surely?
table(H3K27ac_T@seqnames)

##### How much to filter the scores by? #####

length(H3K27ac_T$score[H3K27ac_T$score>0]) #46842450

length(H3K27ac_T$score[H3K27ac_T$score>0.05]) #40190870

length(H3K27ac_T$score[H3K27ac_T$score>mean(H3K27ac_T$score)]) #21670698

hist(H3K27ac_T$score[H3K27ac_T$score>mean(H3K27ac_T$score)])
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


#need to change the seqlengths of the histones to make it work for findOverlapPairs

txdb_seqlengths <- seqlengths(mouse_genes_promoters_coords)

txdb_sl_reordered = txdb_seqlengths[c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")]


normalised_H3K27ac_T <- H3K27ac_T
seqlengths(normalised_H3K27ac_T) <- txdb_sl_reordered


norm_reordered_intersect <- findOverlapPairs(normalised_H3K27ac_T, mouse_genes_promoters_coords)

##### Get the sequences ##### 

H3K27ac_T_promoter_sequences <- getSeq(mouse_genome, norm_reordered_intersect@second)

#histone scores

H3K27ac_T_score <- norm_reordered_intersect@first@elementMetadata$score

hist(H3K27ac_T_score)

max(H3K27ac_T_score) #9.022825

#promoter coordinates

head(norm_reordered_intersect@second@ranges)

promoter_H3K27ac_T_start <- norm_reordered_intersect@second@ranges@start

promoter_H3K27ac_T_end<- norm_reordered_intersect@second@ranges@start + norm_reordered_intersect@second@ranges@width - 1

#histone coordinates

H3K27ac_T_start <- norm_reordered_intersect@first@ranges@start
H3K27ac_T_end <- norm_reordered_intersect@first@ranges@start + norm_reordered_intersect@first@ranges@width - 1

##### Get the context of the histones #####

context_H3K27ac_T_start <- H3K27ac_T_start - promoter_H3K27ac_T_start
context_H3K27ac_T_end <- H3K27ac_T_end - promoter_H3K27ac_T_start

##### IDs for the promoters #####

H3K27ac_T_promoter_chromosome <- norm_reordered_intersect@second@elementMetadata$tx_chrom

H3K27ac_T_promoter_coord_names <- paste0(H3K27ac_T_promoter_chromosome, "_", promoter_H3K27ac_T_start, "_", promoter_H3K27ac_T_end, collapse = ";")

H3K27ac_T_promoter_names <- as.list(unlist(strsplit(H3K27ac_T_promoter_coord_names, "[;]")))

##### Assign each histone to it's promoter region #####

H3K27ac_T_sequences_list <- as.character(H3K27ac_T_promoter_names, use.names=TRUE)
H3K27ac_T_sequence_names <- as_tibble(H3K27ac_T_sequences_list)

H3K27ac_T_tibble <- H3K27ac_T_sequence_names %>% 
  add_column(as.numeric(context_H3K27ac_T_start)) %>% 
  add_column(as.numeric(context_H3K27ac_T_end)) %>% 
  add_column(as.character(H3K27ac_T_promoter_sequences)) %>% 
  add_column(H3K27ac_T_promoter_names) %>% 
  add_column(as.numeric(H3K27ac_T_score)) %>% 
  dplyr::rename(chromosome = value) %>%
  dplyr::rename(H3K27ac_T_start = `as.numeric(context_H3K27ac_T_start)`) %>%
  dplyr::rename(H3K27ac_T_end = `as.numeric(context_H3K27ac_T_end)`) %>%
  dplyr::rename(seq = `as.character(H3K27ac_T_promoter_sequences)`) %>% 
  dplyr::rename(scores = `as.numeric(H3K27ac_T_score)`)

#need to process it somehow so that we get a feature vector and collapse the duplicatesssssss

min(H3K27ac_T_tibble$H3K27ac_T_start) #-49
max(H3K27ac_T_tibble$H3K27ac_T_end) #2798

#convert all the negative starts to 0s and too large ends to 2749

H3K27ac_T_fiddle <- H3K27ac_T_tibble
H3K27ac_T_fiddle$H3K27ac_T_start <- replace(H3K27ac_T_fiddle$H3K27ac_T_start, H3K27ac_T_fiddle$H3K27ac_T_start<0, 0)
H3K27ac_T_fiddle$H3K27ac_T_end <- replace(H3K27ac_T_fiddle$H3K27ac_T_end, H3K27ac_T_fiddle$H3K27ac_T_end > 2749, 2749)  

# get a column with the histone range  
#H3K27ac_T_fiddle <-  H3K27ac_T_fiddle %>% 
  # mutate(histone_range = mapply(function(x,y) seq(from = x, to = y, by = 1), x = .$H3K27ac_T_start, y = .$H3K27ac_T_end ))

# collapse all the histones for each promoter seqeunce
H3K27ac_T_processed <-  H3K27ac_T_fiddle %>%
  filter(H3K27ac_T_score > mean(H3K27ac_T_score))


#FV

H3K27ac_T_bare_minimum <- H3K27ac_T_processed %>% 
  dplyr::select(H3K27ac_T_start, H3K27ac_T_end, scores)

H3K27ac_T_bare_minimum  <- as.data.table(H3K27ac_T_bare_minimum)

H3K27ac_T_fv <- as.data.frame(matrix(nrow = 1, ncol = 2750))
H3K27ac_T_fv <- as.numeric(H3K27ac_T_fv)

H3K27ac_T_replace <- function(i) {
  H3K27ac_T_fv[((H3K27ac_T_bare_minimum$H3K27ac_T_start[i])+1):((H3K27ac_T_bare_minimum$H3K27ac_T_end[i]+1))] <- H3K27ac_T_bare_minimum$scores[i]
  H3K27ac_T_fv
}


dim(H3K27ac_T_bare_minimum) #615894


started.at=proc.time()
H3K27ac_T_temp <- lapply(1:615894, H3K27ac_T_replace)
cat("Finished in",timetaken(started.at),"\n") #28.3 seconds whattt 

peek1 <- H3K27ac_T_temp[1:6]

H3K27ac_T_processed_seq <- H3K27ac_T_processed %>% 
  dplyr::select(seq) %>% 
  group_by(seq)

H3K27ac_T_indices <- group_rows(H3K27ac_T_processed_seq)
H3K27ac_T_coalesced <- vector(length = 1)

start <- lapply(H3K27ac_T_indices, min)
hold <- lapply(H3K27ac_T_indices, length)
end <- mapply(function(x,y)list(sum(x, y, -1)),x = start, y = hold)


H3K27ac_T_coalesce <- function(i) {
  H3K27ac_T_coalesced <- list(as.numeric(c(fcoalesce(H3K27ac_T_temp[c(H3K27ac_T_indices[[i]])]))))
  H3K27ac_T_coalesced
}

length(H3K27ac_T_indices) #19990

started.at=proc.time()
H3K27ac_T_matrix_fv <- sapply(1:19990, H3K27ac_T_coalesce)
cat("Finished in",timetaken(started.at),"\n") #4.952 seconds

H3K27ac_T_matrix_fv <- as.data.table(t(as.data.frame(H3K27ac_T_matrix_fv )))

H3K27ac_T_matrix_fv[is.na(H3K27ac_T_matrix_fv)] <- 0

peek <- H3K27ac_T_matrix_fv[1:6]



#save what I have painstakingly made

H3K27ac_T_sorted <- H3K27ac_T_processed %>% 
  arrange(H3K27ac_T_processed$seq) %>% 
  distinct(seq, .keep_all = T)

H3K27ac_T_sorted <- H3K27ac_T_sorted %>% 
  dplyr::select(chromosome, seq)

H3K27ac_T_final <- H3K27ac_T_matrix_fv %>% 
  add_column(H3K27ac_T_sorted, .before = 1)

save.image("H3K27ac_T.RData")

setwd("/home/cjls4/feature_vectors/")
save(H3K27ac_T_matrix_fv, file = "H3K27ac_T_matrix_fv.RData")

save(H3K27ac_T_final, file = "H3K27ac_T_final")
