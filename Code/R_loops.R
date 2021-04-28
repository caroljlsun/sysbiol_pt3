##### things in R-loops directory#####

#instructions.md
# are simple instructions on where to get the following files

#mm10.rloop.filtered.nameEdit.bed
# has columns
# chr start end mm10.chr.?.?.m? , 0?, + or - (strand), start end, something, something, something*2, something*2

#mm10.rloop.filtered.mergedRegion.bed 
# chrs, start, end, more chromosome info (unique name perhaps), "0", strand

load("")

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

##### Global options #####

setwd(dir = "/storage/CarolSun/R-Loops//")
options(scipen=999)
options(expressions = 500000)

##### Importing data #####

loops <- readGeneric("mm10.rloop.filtered.mergedRegion.bed",
                     chr = 1,
                     start = 2,
                     end = 3,
                     strand = 6)
head(loops)


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

intersect_loops <- findOverlapPairs(loops, mouse_genes_promoters_coords, type = "within")


#get sequences

sequences_loops <- getSeq(mouse_genome, intersect_loops@second)

#get promoter coordinates
head(intersect_loops@second@ranges)

promoter_intersect_loops_start <- intersect_loops@second@ranges@start

promoter_intersect_loops_end<- intersect_loops@second@ranges@start + intersect_loops@second@ranges@width - 1

#get loop coordinates

loop_start <- intersect_loops@first@ranges@start
loop_end <- intersect_loops@first@ranges@start + intersect_loops@first@ranges@width - 1


# make contextual loop coordinates
# so that we get the start and ends of the loops in the context of the promoters themselves
# with the beginning being the start of the promoter etc.

context_loop_start <- loop_start - promoter_intersect_loops_start
context_loop_end <- loop_end - promoter_intersect_loops_start

#make a featurevector

# first get some IDs for the promoters

loop_promoter_chromosome <- intersect_loops@second@elementMetadata$tx_chrom

loop_promoter_coord_names <- paste0(loop_promoter_chromosome, "_", promoter_intersect_loops_start, "_", promoter_intersect_loops_end, collapse = ";")

loop_promoter_names <- as.list(unlist(strsplit(loop_promoter_coord_names, "[;]")))

# they're not unique though, so may need to use tx_ids as well

loop_unique_names <- intersect_loops@second@elementMetadata$tx_name

##### assignment of specific loop positions to each promoter #####

loop_sequences_list <- as.character(loop_promoter_names, use.names=TRUE)
loop_sequence_names <- as_tibble(loop_sequences_list)

loop_tibble <- loop_sequence_names %>% 
  add_column(context_loop_start) %>% 
  add_column(context_loop_end) %>% 
  add_column(as.character(sequences_loops)) %>% 
  add_column(loop_unique_names) %>% 
  dplyr::rename(chromosome = value) %>%
  dplyr::rename(loop_start = context_loop_start) %>%
  dplyr::rename(loop_end = context_loop_end) %>%
  dplyr::rename(seq = `as.character(sequences_loops)`)

##### featurevector #####

loop_fv <- matrix( data = 0, nrow = length(context_loop_start), 
                      ncol = 2750)

for (i in 1:dim(loop_fv)[1]){
  loop_fv[i, (loop_tibble$loop_start[i]:loop_tibble$loop_end[i])] <- 1
}

#check

table(loop_fv[1,])
table(loop_fv[600,])

#####
loop_fv_1 <- loop_fv
loop_fv_1[loop_fv_1 == 0] <- NA
loop_fv_1 <- data.table(loop_fv_1)
loop_fv_2 <- lapply(loop_fv_1, as.numeric)
loop_fv_3 <- data.table::transpose(loop_fv_2)

loops_processed_seq <- loop_tibble %>% 
  dplyr::select(seq) %>% 
  group_by(seq)

loops_indices <- group_rows(loops_processed_seq)
loops_coalesced <- vector(length = 1)

does_it_work <- fcoalesce(loop_fv_3[c(1,2)])

loops_coalesce <- function(i){
  loops_coalesced <- list(as.numeric(c(fcoalesce(loop_fv_3[c(loops_indices[[i]])]))))
  loops_coalesced
}

length(loops_indices) # 8865

loops_matrix_fv <- sapply(1:8865, loops_coalesce)

loops_matrix <- ldply(loops_matrix_fv)

loops_matrix[is.na(loops_matrix)] <- 0

# prepare to save

sorted_loops <- loop_tibble %>% 
  arrange(loop_tibble$seq) %>% 
  distinct(seq, .keep_all = T)

sorted_loops <- sorted_loops %>% 
  dplyr::select(chromosome, seq)

loops_fv_final <- loops_matrix %>% 
  add_column(sorted_loops, .before = 1)

save.image("loops.RData")

setwd("/home/cjls4/feature_vectors/")

save(loops_fv_final, file = "loops_fv_final")

