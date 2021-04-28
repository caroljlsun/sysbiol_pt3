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


#yeoh so something is very wrong with the above method
#can probably adapt the group rows function as well to use here...???


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



BT_seq <- BT_tibble %>% 
  dplyr::select(seq) %>% 
  group_by(seq)

BT_indices <- group_rows(BT_seq)
BT_coalesced <- vector(length = 1)


BT_coalesce <- function(i) {
  BT_coalesced <- list(as.numeric(c(fcoalesce(BT_temp[c(BT_indices[[i]])]))))
  BT_coalesced
}

length(BT_indices) #4867

BT_matrix_fv_corrected<- sapply(1:4867, BT_coalesce)

BT_matrix_fv_corrected_matrix<- as.data.table(t(as.data.frame(BT_matrix_fv_corrected)))

View(BT_matrix_fv_corrected[[7]]) #should have values at 1728 and 1733

BT_matrix_fv_corrected_matrix[is.na(BT_matrix_fv_corrected_matrix)] <- 0



# BT_processed_hmm <-  BT_tibble %>%
#   group_by(seq) %>% 
#   transmute(totalpositions = list(unique(positions+1)),
#             seq = seq) %>%
#   distinct(seq, .keep_all = T)
  
BT_processed_2<-  BT_tibble %>%
  group_by(seq) %>% 
  dplyr::mutate(totalpositions = list(unique(positions+1))) %>%
  distinct(seq, .keep_all = T)

BT_processed_3 <- BT_processed_2 %>% 
  filter(BT_percentages>30)

# 
# sorted_BT <- BT_tibble %>% 
#   arrange(BT_tibble$seq) %>% 
#   distinct(seq, .keep_all = T)
# 
# 
# BT_corrected <- sorted_BT %>% 
#   dplyr::select(-positions) %>%
#   add_column(tibble(matrix(BT_indices))) %>% 
#   dplyr::rename(positions = `matrix(BT_indices)`)
# 
# #filter by 30% 
# BT_corrected_1 <-  BT_corrected %>% 
#   filter(BT_percentages>30)


# The next bit is to tidy up the methylation positions, to remove "c()" essentially

BT_processed_3[6] <- lapply(BT_processed_3[6], gsub, pattern=c("c"), replacement='')
BT_processed_3[6] <- lapply(BT_processed_3[6], gsub, pattern=c("\\("), replacement='')
BT_processed_3[6] <- lapply(BT_processed_3[6], gsub, pattern=c("\\)"), replacement='')

# now it turns out you can't just use the data above, as the strands are not sorted
# so repeat the above, but for positive strands only

# BT_corrected_positive <-  BT_corrected_1 %>% 
#   filter(BT_promoter_strand == "+")
# 
# BT_corrected_negative <- BT_corrected_1 %>% 
#   filter(BT_promoter_strand == "-")


BT_corrected_positive_3 <-  BT_processed_3 %>% 
  filter(BT_promoter_strand == "+")



BT_processed_4<-  BT_tibble %>%
  group_by(seq) %>% 
  dplyr::mutate(totalpositions = list(unique(2750-(positions+1)))) %>%
  distinct(seq, .keep_all = T)

BT_processed_5 <- BT_processed_4 %>% 
  filter(BT_percentages>30)

BT_corrected_negative_3 <- BT_processed_5 %>% 
  filter(BT_promoter_strand == "-")

BT_corrected_negative_3[6] <- lapply(BT_corrected_negative_3[6], gsub, pattern=c("c"), replacement='')
BT_corrected_negative_3[6] <- lapply(BT_corrected_negative_3[6], gsub, pattern=c("\\("), replacement='')
BT_corrected_negative_3[6] <- lapply(BT_corrected_negative_3[6], gsub, pattern=c("\\)"), replacement='')

# then repeat again for the negative strands
# since I changed the promoter definition this will break 
# so change it if needed
# BT_processed_negative <-  group_by(BT_tibble, seq) %>%
#   filter(BT_percentages>50) %>% 
#   transmute(totalpositions = list(unique(5500-(positions+1))), #note the 5500-positions+1. This is very important and took me a while to find out.
#             names = unique(BT_promoter_names),
#             percentages = list(unique(BT_percentages)),
#             strand = unique(BT_promoter_strand)) %>%
#   ungroup() %>% 
#   distinct() %>% 
#   filter(strand == "-")

##### try with only positive strands #####

write.fasta( as.list(BT_corrected_positive_3$seq), BT_corrected_positive_3$BT_promoter_names, file.out = "positive_seq_only.fasta", as.string = T)

modified_BT_corrected_positive <- tibble(BT_corrected_positive_3$totalpositions) %>% 
  add_column(BT_corrected_positive_3$BT_promoter_names)

# double all rows 
new_positive <- modified_BT_corrected_positive[rep(1:nrow(modified_BT_corrected_positive),1,each=2),]

# replace all duplicates with appropriate sequence names
new_positive$`BT_corrected_positive_3$positions`[c(seq(1, dim(new_positive)[1], by=2))] <- paste(">", new_positive$`BT_corrected_positive_3$BT_promoter_names`[c(seq(2, dim(new_positive)[1], by=2))])

new_positive$`BT_corrected_positive_3$positions`<- lapply(new_positive$`BT_corrected_positive_3$positions`, gsub, pattern=c(" "), replacement='')

write(unlist(new_positive$`BT_corrected_positive_3$positions`), file = "positive_positions.fasta", sep = ">")

#position mentioned in positive_positions.fasta for sequence chr1_5585993_5588742 is not valid. Letter 'C' was not found at position 2661

# checking where the C's are
chr1_5585993_5588742 <- gregexpr(pattern ='C', BT_corrected_positive_3$seq[[1]])
chr1_5585993_5588742_t <- tibble(unlist(chr1_5585993_5588742))
#there is a c at 5162!


##### try with only negative strands #####

write.fasta( as.list(BT_corrected_negative_3$seq), BT_corrected_negative_3$BT_promoter_names, file.out = "negative_seq_only.fasta", as.string = T)

modified_BT_processed_negative <- tibble(BT_corrected_negative_3$totalpositions) %>% 
  add_column(BT_corrected_negative_3$BT_promoter_names)

# double all rows 
new_negative <- modified_BT_processed_negative[rep(1:nrow(modified_BT_processed_negative),1,each=2),]

# replace all duplicates with appropriate sequence names
new_negative$`BT_corrected_negative_3$totalpositions`[c(seq(1, dim(new_negative)[1], by=2))] <- paste(">", new_negative$`BT_corrected_negative_3$BT_promoter_names`[c(seq(2, dim(new_negative)[1], by=2))])


new_negative$`BT_corrected_negative_3$totalpositions`<- lapply(new_negative$`BT_corrected_negative_3$totalpositions`, gsub, pattern=c(" "), replacement='')

write(unlist(new_negative$`BT_corrected_negative_3$totalpositions`), file = "negative_positions.fasta", sep = ">")

##### Troubleshooting the negative strand #####

#ERROR: position mentioned in negative_positions.fasta for sequence chr1_5916899_5922398 
#is not valid. Letter 'C' was not found at position 456


# checking where the C's are
gregexpr(pattern ='C', BT_processed_negative$seq[[2]])

# there is a C at 449 and 463 ... so nowhere near 456...

chr1_5916899_5922398 <- gregexpr(pattern ='C', BT_processed_negative$seq[[2]])
chr1_5916899_5922398_t <- tibble(unlist(chr1_5916899_5922398))
nchar(BT_processed_negative$seq[2])
5500-456 = 5044
#indeed there is a c at 5044!

chr1_4785227_4790726 <- gregexpr(pattern = 'C', BT_processed_negative$seq[[1]])
chr1_4785227_4790726_t <- tibble(unlist(chr1_4785227_4790726))
5500-296 = 5204
#indeed there is a c at 5204!

chr1_12990464_12995963 <- gregexpr(pattern = 'C', BT_processed_negative$seq[[3]])
chr1_12990464_12995963_t <- tibble(unlist(chr1_12990464_12995963))
# positions are 237, 256, 306, 348
# there's a c at 306 only
# if we reverse comp them...
5500-237 = 5263
5500-256 = 5244
5500-306 = 5194
5500-348 = 5152

#so the rev positions are 5263, 5244, 5194, 5152
# yes!!! there are c's in all of these reverse comp positions!!!

chr1_4785477_4788226 <- gregexpr(pattern = 'C', BT_corrected_negative_3$seq[[1]])
chr1_4785227_4790726_t <- tibble(unlist(chr1_4785477_4788226))
# should be at position 46

2750-46

##### DNASHAPER #####

pred_methy_pos <- getShape("positive_seq_only.fasta", shapeType = "All",  methylate = TRUE, methylatedPosFile = "positive_positions.fasta")
plotShape(pred_methy_pos$MGW)
plotShape(pred_methy_pos$)


pred_methy_neg <- getShape("negative_seq_only.fasta", methylate = TRUE, methylatedPosFile = "negative_positions.fasta")
plotShape(pred_methy_neg$MGW)

# MGW can be replaced with HelT, Rise, Roll, Shift, Slide, Tilt, Buckle, Opening, ProT, Shear, Stagger, Stretch or EP
plotShape(pred_methy_pos$HelT)
plotShape(pred_methy_pos$ProT)
plotShape(pred_methy_pos$Rise)

#type_final <- c("1-MGW", "1-ProT", "1-Roll", "1-HelT", "1-Rise", "1-Shift", "1-Slide", "1-Tilt", "1-Buckle", "1-Opening", "1-Shear", "1-Stagger", "1-Stretch", "1-EP")

# plotShape(pred_methy$Tilt) ### Error in colMeans(shapeMatrix, na.rm = TRUE) : 
# 'x' must be an array of at least two dimensions
# In addition: There were 15 warnings (use warnings() to see them)

heatShape(pred_methy_pos$ProT, 10)

### will need to change the data shape before doing this but I think it will be a good diagram
# *Note that the input data should only contain one sequence.
# fn2 <- system.file("extdata", "SingleSeqsample.fa", package = "DNAshapeR")
# pred2 <- getShape(fn2)
# trackShape(fn2, pred2) # Only for single sequence files

##### Export the DNAshapes #####

pred_meth_p_chr1 <- getShape("positive_seq_chr1.fasta", methylate = TRUE, methylatedPosFile = "positive_positions_chr1.fasta")
plotShape(pred_meth_p_chr1$MGW)

### Feature vectors/encoding 

type1<- c("1-shape")
vector1 <- encodeSeqShape("positive_seq_only.fasta", pred_methy_pos, type1)

u1 <- data.table(vector1) #1524 10986, 3048/2=1524


type2 <- c("1-MGW")
vector2 <- encodeSeqShape("positive_seq_only.fasta", pred_methy_pos, type2)

u2 <- data.table(vector2) # 2746 variables. It seems like the feature made ignores the first and last 2 bases

# type_final <- c(“n-MGW”, “n-ProT”, “n-Roll”, “n-HelT”, “n-Rise”, “n-Shift”, “n-Slide”,
#                  “n-Tilt”, “n-Buckle”, “n-Opening”, “n-Shear”, “n-Stagger”, “n-Stretch”, “n-EP"")

#type_final <- c("1-MGW", "1-ProT", "1-Roll", "1-HelT", "1-Rise", "1-Shift", "1-Slide", "1-Tilt", "1-Buckle", "1-Opening", "1-Shear", "1-Stagger", "1-Stretch", "1-EP")

# all the positive seqs

vector_pos_MGW <- encodeSeqShape("positive_seq_only.fasta", pred_methy_pos, c("1-MGW"))
vector_final_pos_MGW <- as.data.table(vector_pos_MGW) %>% 
  add_column(BT_corrected_positive_3$seq, .before = 1) %>% 
  dplyr::rename(seq = `BT_corrected_positive_3$seq`)

vector_pos_ProT <- encodeSeqShape("positive_seq_only.fasta", pred_methy_pos, c("1-ProT"))
vector_final_pos_ProT <- as.data.table(vector_pos_ProT) %>% 
  add_column(BT_corrected_positive_3$seq, .before = 1) %>% 
  dplyr::rename(seq = `BT_corrected_positive_3$seq`)

vector_pos_Roll <- encodeSeqShape("positive_seq_only.fasta", pred_methy_pos, c("1-Roll"))
vector_final_pos_Roll <- as.data.table(vector_pos_Roll) %>% 
  add_column(BT_corrected_positive_3$seq, .before = 1) %>% 
  dplyr::rename(seq = `BT_corrected_positive_3$seq`)

vector_pos_HelT <- encodeSeqShape("positive_seq_only.fasta", pred_methy_pos, c("1-HelT"))
vector_final_pos_HelT <- as.data.table(vector_pos_HelT) %>% 
  add_column(BT_corrected_positive_3$seq, .before = 1) %>% 
  dplyr::rename(seq = `BT_corrected_positive_3$seq`)

# all the negative seqs

vector_neg_MGW <- encodeSeqShape("negative_seq_only.fasta", pred_methy_neg, c("1-MGW"))
vector_final_neg_MGW <- as.data.table(vector_neg_MGW) %>% 
  add_column(BT_corrected_negative_3$seq, .before = 1) %>% 
  dplyr::rename(seq = `BT_corrected_negative_3$seq`)

vector_neg_ProT <- encodeSeqShape("negative_seq_only.fasta", pred_methy_neg, c("1-ProT"))
vector_final_neg_ProT <- as.data.table(vector_neg_ProT) %>% 
  add_column(BT_corrected_negative_3$seq, .before = 1) %>% 
  dplyr::rename(seq = `BT_corrected_negative_3$seq`)

vector_neg_Roll <- encodeSeqShape("negative_seq_only.fasta", pred_methy_neg, c("1-Roll"))
vector_final_neg_Roll <- as.data.table(vector_neg_Roll) %>% 
  add_column(BT_corrected_negative_3$seq, .before = 1) %>% 
  dplyr::rename(seq = `BT_corrected_negative_3$seq`)

vector_neg_HelT <- encodeSeqShape("negative_seq_only.fasta", pred_methy_neg, c("1-HelT"))
vector_final_neg_HelT <- as.data.table(vector_neg_HelT) %>% 
  add_column(BT_corrected_negative_3$seq, .before = 1) %>% 
  dplyr::rename(seq = `BT_corrected_negative_3$seq`)



save.image("Methylkit_2.RData")

setwd("/home/cjls4/feature_vectors/")

save(vector_final_pos_MGW, file = "vector_final_pos_MGW")
save(vector_final_pos_ProT, file = "vector_final_pos_ProT")
save(vector_final_pos_Roll, file = "vector_final_pos_Roll")
save(vector_final_pos_HelT, file = "vector_final_pos_HelT")

save(vector_final_neg_MGW, file = "vector_final_neg_MGW")
save(vector_final_neg_ProT, file = "vector_final_neg_ProT")
save(vector_final_neg_Roll, file = "vector_final_neg_Roll")
save(vector_final_neg_HelT, file = "vector_final_neg_HelT")


load("Methylkit_2.RData")

# the dnashapes don't make values for every base, annoyingly. Even worse, it's not consistent.
# the MGW and ProT ones remove 2 bases from each side
# the Roll and HelT something else, not known yet


dim(vector_final_pos_HelT) #1524 2748

dim(vector_final_pos_MGW) # 1524 2747


assess_bases_MGW <- encodeSeqShape("positive_seq_only.fasta", pred_methy_pos, c("1-mer", "1-MGW"))

head(assess_bases_MGW[1,])

dim(assess_bases_MGW) # 1524 13746

# 4*2750 = 11000 which is the hot encoding
# meaning the MGW has been further reduced to 2746

assess_bases_HelT <- encodeSeqShape("positive_seq_only.fasta", pred_methy_pos, c("1-mer", "1-HelT"))

dim(assess_bases_HelT)# 1524 13747

head(assess_bases_HelT[,13747])

head(vector_final_pos_HelT[,2748])

head(assess_bases_MGW[,13746])

head(vector_final_pos_MGW[,2747])

# how to get the next level shapes?

f1 <- encodeNstOrderShape(1, pred_methy_pos, "EP")

# I don't think I can lmao

##### Stealing DNAshaper's code since the library doesn't want to work #####

encodeNstOrderShape <- function( n, shapeMatrix, shapeType ){
  # assign average value to NA
  #shapeMatrix[ is.na( shapeMatrix ) ] <- 0
  switch( shapeType,
          MGW = { shapeMatrix[ is.na( shapeMatrix ) ] <- 5.072 },
          ProT = { shapeMatrix[ is.na( shapeMatrix ) ] <- -6.792 },
          Roll = { shapeMatrix[ is.na( shapeMatrix ) ] <- -0.698 },
          HelT = { shapeMatrix[ is.na( shapeMatrix ) ] <- 34.326 },
          EP = { shapeMatrix[ is.na( shapeMatrix ) ] <- -6.505 },
          
          Stretch = { shapeMatrix[ is.na( shapeMatrix ) ] <- -0.028 },
          Tilt = { shapeMatrix[ is.na( shapeMatrix ) ] <- -0.008 },
          Buckle = { shapeMatrix[ is.na( shapeMatrix ) ] <- 0.097 },
          Shear = { shapeMatrix[ is.na( shapeMatrix ) ] <- -0.009 },
          Opening = { shapeMatrix[ is.na( shapeMatrix ) ] <- -0.24 },
          Rise = { shapeMatrix[ is.na( shapeMatrix ) ] <- 3.342 },
          Shift = { shapeMatrix[ is.na( shapeMatrix ) ] <- 0.008 },
          Stagger = { shapeMatrix[ is.na( shapeMatrix ) ] <- -0.013 },
          Slide = { shapeMatrix[ is.na( shapeMatrix ) ] <- -1.526 }
  )
  
  
  singleSeq <- FALSE
  if( nrow(shapeMatrix)[1] == 1 )
    singleSeq <- TRUE
  
  # trim both 2 bps end
  if( shapeType == "MGW" || shapeType == "ProT" || shapeType == "EP" ||
      shapeType == "Stretch" || shapeType == "Buckle" ||
      shapeType == "Shear" || shapeType == "Opening" ||
      shapeType == "Stagger" ){
    shapeMatrix <- shapeMatrix[, -c(1, 2, ncol( shapeMatrix )-1,
                                    ncol( shapeMatrix ))]
    
  }else if( shapeType == "Roll" || shapeType == "HelT" ||
            shapeType == "Tilt" || shapeType == "Rise" ||
            shapeType == "Shift" || shapeType == "Slide" ){
    shapeMatrix <- shapeMatrix[, -c(1, ncol( shapeMatrix ))]
  }
  
  if(singleSeq)
    shapeMatrix <- t(shapeMatrix)
  
  # encode k-st feature
  featureVector <- NULL
  if( n == 1 ){
    featureVector = shapeMatrix
    
  }else{
    m <- ncol( shapeMatrix )
    # normalization
    shapeMatrix <- normalizeShape( featureVector = shapeMatrix,
                                   thOrder = 1, shapeType = shapeType,
                                   normalize = TRUE )
    
    for ( i in 1 : ( m-n+1 )){
      feature <- shapeMatrix[, i]
      for ( j in ( i+1 ) : ( i+n-1 ) )
        feature <- feature * shapeMatrix[, j]
      
      featureVector <- cbind( featureVector, unlist( feature ) )
    }
  }
  
  return ( featureVector )
}



##### DNAshaper using data split into chromosomes #####

#using just the positive strand data

BT_processed_positive_chr1 <- BT_processed_positive %>% 
  filter(str_detect(BT_processed_positive$names, "chr1_"))

write.fasta( as.list(BT_processed_positive_chr1$seq), BT_processed_positive_chr1$names, file.out = "positive_seq_chr1.fasta", as.string = T)

modified_BT_processed_positive_chr1 <- tibble(BT_processed_positive_chr1$totalpositions) %>% 
  add_column(BT_processed_positive_chr1$names)

# double all rows 
new_positive <- modified_BT_processed_positive_chr1[rep(1:nrow(modified_BT_processed_positive_chr1),1,each=2),]

# replace all duplicates with appropriate sequence names
new_positive$`BT_processed_positive_chr1$totalpositions`[c(seq(1, dim(new_positive)[1], by=2))] <- paste(">", new_positive$`BT_processed_positive_chr1$names`[c(seq(2, dim(new)[1], by=2))])
new_positive$`BT_processed_positive_chr1$totalpositions`<- lapply(new_positive$`BT_processed_positive_chr1$totalpositions`, gsub, pattern=c(" "), replacement='')

write(unlist(new_positive$`BT_processed_positive_chr1$totalpositions`), file = "positive_positions_chr1.fasta", sep = ">")

pred_meth_p_chr1 <- getShape("positive_seq_chr1.fasta", methylate = TRUE, methylatedPosFile = "positive_positions_chr1.fasta")
plotShape(pred_meth_p_chr1$MGW)

### Feature vectors/encoding 

featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape("positive_seq_chr1.fasta", pred_meth_p_chr1, featureType)


# uhhhh try again, but no more column collapsing
set.seed(12345)

x_bis <- data.frame(featureVector)

y <- rbinom(1445,1, 0.5)
y <- tibble(y)

df_bis <- tibble(y, x_bis)

splitter <- sample(c(rep(0, round(0.67*nrow(df_bis), 0)),
                     rep(1, round(0.33*nrow(df_bis), 0))))
splitter
table(splitter)

df_bis_train <- df_bis[splitter == 0, ]
df_bis_test <- df_bis[splitter == 1, ]

xbistrain <-as.matrix(df_bis_train[2:13745])
ybistrain <-as.matrix(df_bis_train[1])

xbistest <- as.matrix(df_bis_test[2:13745])
ybistest <- as.matrix(df_bis_test[1])

#dimensions declarations

dim(x_bis)
#1445x13746

#keras troubleshooting

reticulate::py_available()
reticulate::import("keras.models")

reticulate::py_config()

# Initialize a sequential model
mod <- Sequential() 

# Add layers to the model
mod$add(Flatten(input_shape = 13744))
mod$add(Activation("relu"))
mod$add(Dense(100))
mod$add(Activation("relu"))
mod$add(Dense(1))
mod$add(Activation("sigmoid"))
mod

keras_compile(mod,  loss = 'binary_crossentropy', metrics = c('binary_accuracy'), optimizer = 'SGD')

set.seed(12345)

keras_fit(mod, xbistrain, ybistrain, validation_split = 0.1, epochs=15, batch_size=100, verbose = 1)

deepviz::plot_model(mod)

probY <- keras_predict_proba(mod, xbistest)
predictY <- keras_predict_classes(mod, xbistest)

evaluate(mod, xbistest, ybistest)
#x can be a matrix so no need for the  I did