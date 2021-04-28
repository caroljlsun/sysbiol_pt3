#load workspace

#load("Expression_workspace.RData")

# libraries

library(limma)
library(data.table)
library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(dplyr)
library(GenomicFeatures)
library(genomation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(useful)
library(biomaRt)

#Global options
setwd(dir = "/storage/CarolSun/Expression")
options(scipen=999)
options(expressions = 500000)
#browseURL('https://www.youtube.com/watch?v=QH2-TGUlwu4')

#CPM = counts per million
CPM <- read_tsv("Blueprint_Expression_CPM_annotated.tsv")
CPM_ranges <- CPM %>% 
  dplyr::select(chromosome_name:strand, ensembl_gene_id:external_gene_name) %>% 
  rename(start = start_position,
         end = end_position)

CPM_ranges$chromosome_name <- paste0("chr", CPM_ranges$chromosome_name)

CPM_ranges$strand[CPM_ranges$strand==1] <- as.character("+")
CPM_ranges$strand[CPM_ranges$strand==-1] <- as.character("-")

#CPM_metadata <- CPM %>% 
 # dplyr::select(ensembl_gene_id:tcell_M)
#The data CPM has ensembl gene IDs. The txdb stuff uses entrez ids
#Here I tried to use biomart to get the CPM data to use entrez ids, then find the intersect...


CPM_keys <- CPM %>% 
  dplyr::select(ensembl_gene_id)

ensembl <- useMart("ensembl")

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

genes <- CPM_keys$ensembl_gene_id

IDs_CPM <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "entrezgene_id"), values = c(genes), mart = mart)

#turn the genes positions into a granges object

CPM_ranges <- as_granges(CPM_ranges, seqnames = chromosome_name)

#then intersect with promoter regions???

mouse_annotation <- TxDb.Mmusculus.UCSC.mm10.knownGene
mouse_genome <- BSgenome.Mmusculus.UCSC.mm10

seqlevels(mouse_annotation) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                 "chr11","chr12","chr13","chr14","chr15", "chr16","chr17","chr18","chr19","chrX","chrY")

geneID_mouse_genes <- keys(mouse_annotation, keytype = "GENEID")

#find the shared subset of entrez gene ids

geneIDs_CPM <- tibble(IDs_CPM) %>% 
  dplyr::select(entrezgene_id) %>% 
  rename(gene_id = entrezgene_id)


geneIDs_mouse_genes <- as.numeric(as.character(geneID_mouse_genes)) %>% 
  tibble()%>% 
  rename(gene_id = ".")

#why is this broken ffssssssssssss

intersect_geneIDs <-  dplyr::intersect(geneIDs_CPM, geneIDs_mouse_genes)

#get the intersection seq ranges

#need to get the corresponding ensembl ids from the gene ids

#first, the matching gene ids from the intersection and the IDs_CPM list

intersect_ensembl_IDs <-  inner_join(IDs_CPM, intersect_geneIDs, by = c("entrezgene_id" = "gene_id"))

intersect_IDs_ranges <- CPM_ranges[(elementMetadata(CPM_ranges)[,1] %in% intersect_ensembl_IDs$ensembl_gene_id)]


#get the standard promoter regions

mouse_genes_transcript_coords <- transcriptsBy(mouse_annotation, by = "gene")

mouse_genes_promoters <- getPromoterSeq(mouse_genes_transcript_coords, mouse_genome, upstream = 2500, downstream = 250)

mouse_genes_promoters_coords   <- promoters(mouse_annotation, upstream=2500, downstream=250, 
                                            columns=c("tx_id", "tx_chrom", "tx_strand", "tx_name", "gene_id"))

intersect_promoters_IDs <- findOverlapPairs(intersect_IDs_ranges, mouse_genes_promoters_coords)

#Choice, do I do the standard or do I do one value only per promoter?
#He said to just do one value per promoter...

head(intersect_promoters_IDs@first)

#get sequences

sequences_intersect <- getSeq(mouse_genome, intersect_promoters_IDs@second)
head(sequences_intersect)

#get promoter coordinates
head(intersect_promoters_IDs@second@ranges)

promoter_intersect_start <- intersect_promoters_IDs@second@ranges@start
promoter_intersect_end <- intersect_promoters_IDs@second@ranges@start + intersect_promoters_IDs@second@ranges@width - 1

# check width

intersect_promoter_widths <- promoter_intersect_end - promoter_intersect_start + 1 #should this be +1...

intersect_promoter_chromosome <- intersect_promoters_IDs@second@elementMetadata$tx_chrom

intersect_promoter_coord_names <- paste0(intersect_promoter_chromosome, "_", promoter_intersect_start, "_", promoter_intersect_end, collapse = ";")

intersect_promoter_names <- as.list(unlist(strsplit(intersect_promoter_coord_names, "[;]")))


intersect_sequences_list <- as.character(intersect_promoter_names, use.names=TRUE)
intersect_sequence_names <- as_tibble(intersect_sequences_list)

# Get the expression values for the B and T cells

bcellf <- intersect_promoters_IDs@first@elementMetadata$bcell_F

bcellm <- intersect_promoters_IDs@first@elementMetadata$bcell_M

tcellf <- intersect_promoters_IDs@first@elementMetadata$tcell_F

tcellm <- intersect_promoters_IDs@first@elementMetadata$tcell_M

# get the gene names too

gene_names <- intersect_promoters_IDs@first@elementMetadata$external_gene_name

#Collapse into a tibble or something

intersect_tibble <-intersect_sequence_names %>%
  add_column(as.character(sequences_intersect)) %>% 
  add_column(bcellf) %>% 
  add_column(bcellm) %>% 
  add_column(tcellf) %>% 
  add_column(tcellm) %>% 
  add_column(gene_names) %>% 
  rename(names = value) %>% 
  rename(seq = `as.character(sequences_intersect)`)%>%
  distinct(., seq, .keep_all = TRUE) %>% 
  dplyr::filter(intersect_tibble$bcellf != 0 & intersect_tibble$bcellm != 0 & intersect_tibble$tcellf != 0 & intersect_tibble$tcellm != 0)

#should be 23264

table(intersect_tibble$bcellf == 0 & intersect_tibble$bcellm == 0 & intersect_tibble$tcellf == 0 & intersect_tibble$tcellm == 0)

#save

save.image("Expression_workspace.RData")


setwd("/home/cjls4/feature_vectors/")

save(intersect_tibble, file = "expression_tibble.RData")
