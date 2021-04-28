load()

##### Loading libraries #####

library(methylKit)
library(GenomicFeatures)
library(GenomicRanges)
library(genomation)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DNAshapeR)
library(seqinr)
library(fields)
library(Biostrings)
library(useful)


##### Global options to make sure the code actually works #####
setwd(dir = "/storage/CarolSun/Methylation/")
options(scipen=999)
options(expressions = 500000)

##### defining functions #####

### read_bismark_processed
# This function will convert the methylation files we currently have into a
# methyKit object, specifically a methylRawList
# which can then be manipulated into other forms

read_bismark_processed<-function( location,sample.id,assembly="unknown",treatment,
                                  context="CpG",min.cov=10)
{
  if(length(location)>1){
    stopifnot(length(location)==length(sample.id),
              length(location)==length(treatment))
  }
  
  result=list()
  for(i in 1:length(location)){
    df=read_tsv(location[[i]])
    
    
    # make the object (arrange columns of df), put it in a list
    result[[i]]= new("methylRaw",data.frame(chr=df[,1],start=df[,2],end=df[,2],
                                            strand="*",coverage=(df[,3]+df[,4]),
                                            numCs=df[,3],numTs=df[,4]),
                     sample.id=sample.id[[i]],
                     assembly=assembly,context=context,resolution="base"
    )
    
  }
  
  if(length(result) == 1){
    return(result[[1]])
  }else{
    
    new("methylRawList",result,treatment=treatment)
    
  }
  
}
# argument: location can be a list of file names
# argument: treatment can be a list of binary numbers, c(0,0,1,1)
# argument: sample.id should be a list of characters, list("a", "b", "c")

### rename_methyl
# This function will rename the names of the columns of the methylRawList
# You need to do this to change this object into methylRawDB's
rename_methyl <- function(object)
{
  propernames <- list("chr", "start", "end", "strand", "coverage", "numCs", "numTs")
  
  for(i in 1:length(propernames)){
    object@names[i] <- propernames[[i]]
  }
  return(object)
}
# argument: object is your methylRawList.
# propernames: if needed, can be changed.

#####  (Ignore)Importing sample data - OX only #####

# let's get several .tsv's converted
# and use only the "OX" types!!!!
# this will be a small group, which acts as a quality control/ check
# thus the name, qc_group

sample_ids <- list("sample_F_B_9", "sample_F_B_11", "sample_F_T_14", "sample_F_T_16")
treatments <- c(1,1,0,0) # where 1 = B cell, 0 = T cell

qc_group_tsv <- list("B6_F_B_OX.9.methy.combined_strand.mincov_10.tsv", 
                     "B6_F_B_OX.11.methy.combined_strand.mincov_10.tsv",
                     "B6_F_T_OX.14.methy.combined_strand.mincov_10.tsv",
                     "B6_F_T_OX.16.methy.combined_strand.mincov_10.tsv")

##### (Ignore)Wrangle #####

qc_group_methyl <- read_bismark_processed(qc_group_tsv,
                                          sample.id = sample_ids,
                                          treatment = treatments)
#check that it is the right object class
#class(qc_group_methyl)

### renaming the header in all of the samples

for (i in 1:length(qc_group_tsv)){
  qc_group_methyl[[i]] <- rename_methyl(qc_group_methyl[[i]])
}

# check that it worked
# qc_group_methyl[[4]]@names

qc_group_methyl <- as(qc_group_methyl, "methylRawList")

# to make this unite function actually work, you need to make sure that the location numbers 
# are written in expanded form, aka not scientific notation

qc_group_unite <- unite(qc_group_methyl, save.db = T)

# quick check of the data

getMethylationStats(qc_group_methyl[[4]], plot = T)

getCorrelation(qc_group_unite, plot = T)

clusterSamples(qc_group_unite, dist = "correlation", method = "ward", plot = T)

PCASamples(qc_group_unite)

##### (!!!)Import/compile all OX B and T cells #####

samples_BT <- list("F_B_9", "F_B_11", "M_B_2", "M_B_3", "F_T_14", "F_T_16", "M_T_6", "M_T_8")

treatments_BT <- c(1,1,1,1,0,0,0,0)

BT_tsv <- list("B6_F_B_OX.9.methy.combined_strand.mincov_10.tsv",
               "B6_F_B_OX.11.methy.combined_strand.mincov_10.tsv",
               "B6_M_B_OX.2.methy.combined_strand.mincov_10.tsv",
               "B6_M_B_OX.3.methy.combined_strand.mincov_10.tsv",
               "B6_F_T_OX.14.methy.combined_strand.mincov_10.tsv",
               "B6_F_T_OX.16.methy.combined_strand.mincov_10.tsv",
               "B6_M_T_OX.6.methy.combined_strand.mincov_10.tsv",
               "B6_M_T_OX.8.methy.combined_strand.mincov_10.tsv")

# import
BT_methyl <- read_bismark_processed(BT_tsv,
                                    sample.id = samples_BT,
                                    treatment = treatments_BT)
#rename
for (i in 1:length(BT_tsv)){
  BT_methyl[[i]] <- rename_methyl(BT_methyl[[i]])
}

#make db
BT_unite <- methylKit::unite(BT_methyl, save.db = TRUE)

#getMethylationStats(BT_methyl[[2]], plot = T)

#getCorrelation(BT_unite, plot = T)

clusterSamples(BT_unite, dist = "correlation", method = "ward", plot = T)

#PCASamples(BT_unite)

##### (Ignore)Import/compile all OX FB, FT, MB, MT #####



samples_FM <- list("F_B_9", "F_B_11", "F_T_14", "F_T_16", "M_B_2", "M_B_3", "M_T_6", "M_T_8")

treatments_FM <- c(1,1,1,1,0,0,0,0)

FM_tsv <- list("B6_F_B_OX.9.methy.combined_strand.mincov_10.tsv",
               "B6_F_B_OX.11.methy.combined_strand.mincov_10.tsv",
               "B6_F_T_OX.14.methy.combined_strand.mincov_10.tsv", 
               "B6_F_T_OX.16.methy.combined_strand.mincov_10.tsv",
               "B6_M_B_OX.2.methy.combined_strand.mincov_10.tsv",
               "B6_M_B_OX.3.methy.combined_strand.mincov_10.tsv",
               "B6_M_T_OX.6.methy.combined_strand.mincov_10.tsv",
               "B6_M_T_OX.8.methy.combined_strand.mincov_10.tsv")

FM_methyl <- read_bismark_processed(FM_tsv,
                                    sample.id = samples_FM,
                                    treatment = treatments_FM)

# rename
for (i in 1:length(FM_tsv)){
  FM_methyl[[i]] <- rename_methyl(FM_methyl[[i]])
}

#make DB
FM_unite <- methylKit::unite(FM_methyl, save.db = T)

#check out the stats
getMethylationStats(FM_methyl[[6]], plot = T)

getCorrelation(FM_unite, plot = T)

clusterSamples(FM_unite, dist = "correlation", method = "ward", plot = T)

PCASamples(FM_unite)


##### (Quality check, ignore) correlate with genomic positions (specifically the promoter regions) #####

mouse_annotation <- TxDb.Mmusculus.UCSC.mm10.knownGene
mouse_genome <- BSgenome.Mmusculus.UCSC.mm10

#right, going to rip from Russell's code here

my_genes                 <- keys(mouse_annotation)[1:100] # a small selection so far
my_genes_transcript_coords <- transcriptsBy(mouse_annotation, by = "gene")[as.character(my_genes)]

my_genes_promoters     <- getPromoterSeq(my_genes_transcript_coords, mouse_genome, upstream = 2500, downstream = 250)
my_genes_promoters_coords   <- promoters(mouse_annotation, upstream=2500, downstream=250, 
                                         columns=c("tx_id", "tx_chrom", "tx_strand", "tx_name", "gene_id")) 

head(my_genes_promoters)
head(my_genes_promoters_coords)


str(my_genes_promoters_coords@seqinfo)


#test out a smaller set of methylBaseDB
qc_granged <- as(qc_group_unite, "GRanges")

qc_promoters_100_pairs <- findOverlapPairs(qc_granged, my_genes_promoters_coords)
?pintersect()

head(qc_promoters_100_pairs@first)
head(qc_promoters_100_pairs@second)

qc_overlap_only <- qc_promoters_100_pairs@second

trimmed <- trim(qc_overlap_only)
qc_sequences <- getSeq(mouse_genome, trimmed)

### useless code but I dont want to delete it yet ###
# relevant_qc_bool <- str_detect(qc_overlap_only@seqinfo@seqnames, "random", negate = T)
# relevant_qc_names <- c(qc_overlap_only@seqinfo@seqnames[1:21])
# 
# relevant_qc <- qc_overlap_only[seqnames(qc_overlap_only) == "chr2"]
# relevant_chromosomes <- GRanges(seqnames=NULL, ranges=NULL, strand=NULL, seqlengths=NULL, seqinfo=NULL)
# what <- subsetByOverlaps(qc_overlap_only, relevant_chromosomes)




##### (!!!) Use the larger datasets #####

# Generate the promoter sequences of all the mouse genome 

mouse_annotation <- TxDb.Mmusculus.UCSC.mm10.knownGene
mouse_genome <- BSgenome.Mmusculus.UCSC.mm10

seqlevels(mouse_annotation) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                 "chr11","chr12","chr13","chr14","chr15", "chr16","chr17","chr18","chr19","chrX","chrY")

mouse_genes <- keys(mouse_annotation)

mouse_genes_transcript_coords <- transcriptsBy(mouse_annotation, by = "gene")

mouse_genes_promoters <- getPromoterSeq(mouse_genes_transcript_coords, mouse_genome, upstream = 2500, downstream = 250)

mouse_genes_promoters_coords   <- promoters(mouse_annotation, upstream=2500, downstream=250, 
                                            columns=c("tx_id", "tx_chrom", "tx_strand", "tx_name", "gene_id"))

# # check it all looks good
# head(mouse_genes_promoters_coords)
# str(mouse_genes_promoters_coords@seqinfo)
# length(mouse_genes_promoters_coords) #63562
# length(unique(mouse_genes_promoters_coords@ranges)) #45197
# length(mouse_genes_promoters_coords@ranges) #63562

BT_granged <- as(BT_unite, "GRanges")

# length(BT_granged)
# length(unique(BT_granged@ranges))
# head(BT_granged@ranges)

##### (!!!)Find the overlaps between the methylation positions and the promoter coordinates #####

BT_promoters_100_pairs <- findOverlapPairs(BT_granged, mouse_genes_promoters_coords)

# head(BT_promoters_100_pairs@first)
# head(BT_promoters_100_pairs@second)

BT_percentages <- (BT_promoters_100_pairs@first@elementMetadata[[2]]/BT_promoters_100_pairs@first@elementMetadata[[1]])*100

BT_promoter_coords_start <- BT_promoters_100_pairs@second@ranges@start
BT_promoter_coords_end <- BT_promoter_coords_start+BT_promoters_100_pairs@second@ranges@width-1

#get IDs
BT_promoter_chromosome <- BT_promoters_100_pairs@second@elementMetadata$tx_chrom

BT_promoter_coord_names <- paste0(BT_promoter_chromosome, "_" , BT_promoter_coords_start, "_" ,BT_promoter_coords_end, collapse = ";")
BT_promoter_names <- as.list(unlist(strsplit(BT_promoter_coord_names, "[;]")))

#get strands
BT_promoter_strand <- as.list(BT_promoters_100_pairs@second@strand)

#get sequences
BT_overlap_only <- BT_promoters_100_pairs@second

trimmed_BT <- trim(BT_overlap_only)
BT_sequences <- getSeq(mouse_genome, trimmed_BT)


##### (!!!)Calculating methyl postions #####

#first I'll need to get the specific positions of the methylated C in the promoter regions

BT_start <- BT_promoters_100_pairs@second@ranges@start

BT_c_position <- BT_promoters_100_pairs@first@ranges@start

BT_positions <- BT_c_position-BT_start

# assign these specific positions to each promoter region

# extract info from the xtringset

BT_sequences_list <- as.character(BT_sequences, use.names=TRUE)
BT_sequences_names <- attributes(BT_sequences_list) %>% 
  dplyr::as_tibble()

##### (!!!) merge the sequences and positions? into a dataframe #####

BT_tibble <- tibble(BT_sequences_list) %>% 
  add_column(BT_positions) %>% 
  add_column(BT_percentages) %>% 
  add_column(BT_promoter_names) %>% 
  add_column(BT_promoter_strand) %>% 
  dplyr::rename(seq = BT_sequences_list) %>% 
  dplyr::rename(positions = BT_positions) %>% 
  filter(BT_percentages > 0)


#let me check that unique names == unique positions etc

length(unique(BT_tibble$seq)) #4867
length(unique(BT_tibble$positions)) #2682

#I believe that the processes tibble should be 4867 long, then we can "collapse" the positions for duplicated promoters


#also make sure everything is 2750, not 2749

# collapse all the histones for each promoter seqeunce

#FV

BT_try_out <- as.data.frame(BT_tibble)

BT_fv <- as.data.frame(matrix(nrow = 1, ncol = 2750))
BT_fv <- as.numeric(BT_fv)

BT_tibble$positions[1]
BT_tibble$BT_percentages[8]

BT_replace <- function(i) {
  BT_fv[BT_tibble$positions[i]] <- BT_tibble$BT_percentages[i]
  BT_fv
}

dim(BT_tibble)[1]

started.at=proc.time()
BT_temp <- lapply(1:16523, BT_replace)
cat("Finished in",timetaken(started.at),"\n") 

peek1 <- BT_temp[1:6]

#wait it surely needs to be sorted by chr names, before I can do the fun coalesce stuff

BT_seq <- BT_tibble %>% 
  dplyr::select(seq) %>% 
  group_by(seq)

BT_indices <- group_rows(BT_seq)

#check that the indices refer to the right rows
BT_indices[[7]]
look2 <- BT_tibble[13931:13948,]
#It seems like they do

BT_coalesced <- vector(length = 1)

#start <- lapply(BT_indices, min)
#hold <- lapply(BT_indices, length)


huhhhh <- BT_indices[[7]][2] #use this lol?
#I have fucked it with the hold function
# I'm pretty sure this means that it's wrong for all other scripts too fuccccccckkkkkk
#should be easy to fix, but it is a pricy mistake lmao

#no, change the coalesce function so that it coalesces at xyz specific rows, not start:end

please_god <- fcoalesce(BT_temp[c(13931,13948)])
please_god_1 <- fcoalesce(BT_temp[c(1,2,3,4,5,6,7)])
please_god_2 <- fcoalesce(BT_temp[c(BT_indices[[7]])])

#end <- mapply(function(x,y)list(sum(x, y, -1)),x = start, y = hold)

BT_coalesce <- function(i) {
  BT_coalesced <- list(as.numeric(c(fcoalesce(BT_temp[c(BT_indices[[i]])]))))
  BT_coalesced
}

length(BT_indices) #4867

started.at=proc.time()
BT_matrix_fv_corrected<- sapply(1:4867, BT_coalesce)
cat("Finished in",timetaken(started.at),"\n") #0.76 seconds

BT_matrix_fv_corrected_matrix<- as.data.table(t(as.data.frame(BT_matrix_fv_corrected)))

View(BT_matrix_fv_corrected[[7]]) #should have values at 1728 and 1733

BT_matrix_fv_corrected_matrix[is.na(BT_matrix_fv_corrected_matrix)] <- 0

# prepare to save

sorted_methyl <- BT_tibble %>% 
  arrange(BT_tibble$seq) %>% 
  distinct(seq, .keep_all = T)

sorted_methyl <- sorted_methyl %>% 
  dplyr::select(seq)

methyl_fv_final <- BT_matrix_fv_corrected_matrix %>% 
  add_column(sorted_methyl, .before = 1)



#browseURL('https://www.youtube.com/watch?v=QH2-TGUlwu4')

save.image("methylated_positions.RData")


setwd("/home/cjls4/feature_vectors/")
save(BT_matrix_fv_corrected_matrix, file = "BT_methylated_positions.RData")
save(methyl_fv_final, file = "methyl_fv_final")

#need a tibble with the promoter names and sequences please, save that instead of the pure feature matrix