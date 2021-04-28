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
library(caret)
library(keras)
# library(kerasR)
# library(reticulate)
library(tensorflow)
library(useful)


# library(reticulate)
# use_python("/home/cjls4/miniconda3/envs/r-reticulate/bin/python3")
# py_config()
# 
# library(reticulate)
# use_python(python = Sys.which("python3"), required = TRUE)
# py_config()
# 
# use_condaenv("tf")
# tensorflow::install_tensorflow("auto")
# reticulate::py_discover_config()
# 
# tf$constant("Hellow Tensorflow")


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
  rename(seq = BT_sequences_list) %>% 
  rename(positions = BT_positions) %>% 
  filter(BT_percentages > 0)


#let me check that unique names == unique positions etc

length(unique(BT_tibble$seq)) #4867
length(unique(BT_tibble$positions)) #2682

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# I believe that the processes tibble should be 4867 long, then we can "collapse" the positions for duplicated promoters
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# process the dataframe so that you only have sequences with percentage methylations over 50%,
# with unique names, percentages, strands and collapsed methylation positions

# BT_processed <-  group_by(BT_tibble, seq) %>%
#   filter(BT_percentages>50) %>% 
#   transmute(totalpositions = list(unique(positions)),
#             names = unique(BT_promoter_names),
#             percentages = list(unique(BT_percentages)),
#             strand = unique(BT_promoter_strand)) %>%
#   ungroup() %>% 
#   distinct()


BT_corrected

#also make sure everything is 2750, not 2749



# The next bit is to tidy up the methylation positions, to remove "c()" essentially

BT_processed[2] <- lapply(BT_processed[2], gsub, pattern=c("c"), replacement='')
BT_processed[2] <- lapply(BT_processed[2], gsub, pattern=c("\\("), replacement='')
BT_processed[2] <- lapply(BT_processed[2], gsub, pattern=c("\\)"), replacement='')

# now it turns out you can't just use the data above, as the strands are not sorted
# so repeat the above, but for positive strands only

BT_processed_positive <-  group_by(BT_tibble, seq) %>%
  filter(BT_percentages>50) %>% 
  transmute(totalpositions = list(unique(positions+1)), #note the +1. This is very important.
            names = unique(BT_promoter_names),
            percentages = list(unique(BT_percentages)),
            strand = unique(BT_promoter_strand)) %>%
  ungroup() %>% 
  distinct() %>% 
  filter(strand == "+")

BT_processed_positive[2] <- lapply(BT_processed_positive[2], gsub, pattern=c("c"), replacement='')
BT_processed_positive[2] <- lapply(BT_processed_positive[2], gsub, pattern=c("\\("), replacement='')
BT_processed_positive[2] <- lapply(BT_processed_positive[2], gsub, pattern=c("\\)"), replacement='')

# then repeat again for the negative strands
# since I changed the promoter definition this will break 
# so change it if needed
BT_processed_negative <-  group_by(BT_tibble, seq) %>%
  filter(BT_percentages>50) %>% 
  transmute(totalpositions = list(unique(5500-(positions+1))), #note the 5500-positions+1. This is very important and took me a while to find out.
            names = unique(BT_promoter_names),
            percentages = list(unique(BT_percentages)),
            strand = unique(BT_promoter_strand)) %>%
  ungroup() %>% 
  distinct() %>% 
  filter(strand == "-")

BT_processed_negative[2] <- lapply(BT_processed_negative[2], gsub, pattern=c("c"), replacement='')
BT_processed_negative[2] <- lapply(BT_processed_negative[2], gsub, pattern=c("\\("), replacement='')
BT_processed_negative[2] <- lapply(BT_processed_negative[2], gsub, pattern=c("\\)"), replacement='')

##### testout conversion methods to get data into DNAshaper #####

write.fasta( as.list(BT_processed$seq), BT_processed$names, file.out = "test_3.fasta", as.string = T)

modified_BT_processed <- tibble(BT_processed$totalpositions) %>% 
  add_column(BT_processed$names)

# double all rows 
new <- modified_BT_processed[rep(1:nrow(modified_BT_processed),1,each=2),]

# replace all duplicates with appropriate sequence names
new$`BT_processed$totalpositions`[c(seq(1, dim(new)[1], by=2))] <- paste(">", new$`BT_processed$names`[c(seq(2, dim(new)[1], by=2))])
new$`BT_processed$totalpositions`<- lapply(new$`BT_processed$totalpositions`, gsub, pattern=c(" "), replacement='')

#what <- c(seq(1, dim(new)[1], by=2))

length(new$`BT_processed$totalpositions`)
length(unique(new$`BT_processed$totalpositions`))
length(BT_processed$seq)
length(unique(BT_processed$seq))


write.fasta(as.character(paste(BT_processed$totalpositions)), BT_processed$names,
            file.out = "test_3_positions.fasta", as.string = F)

write(unlist(new$`BT_processed$totalpositions`), file = "whyyy.fasta", sep = ">")

positions_fasta <- read.fasta(file = "whyyy.fasta")

#test out in dnashaper

pred_methy <- getShape("test_3.fasta", methylate = TRUE, methylatedPosFile = "whyyy.fasta")
pred_methy$MGW

#nope DNAshaper isnt happy

##### retry with only positive strands #####

write.fasta( as.list(BT_processed_positive$seq), BT_processed_positive$names, file.out = "positive_seq_only.fasta", as.string = T)

modified_BT_processed_positive <- tibble(BT_processed_positive$totalpositions) %>% 
  add_column(BT_processed_positive$names)

# double all rows 
new_positive <- modified_BT_processed_positive[rep(1:nrow(modified_BT_processed_positive),1,each=2),]

# replace all duplicates with appropriate sequence names
new_positive$`BT_processed_positive$totalpositions`[c(seq(1, dim(new_positive)[1], by=2))] <- paste(">", new_positive$`BT_processed_positive$names`[c(seq(2, dim(new_positive)[1], by=2))])
new_positive$`BT_processed_positive$totalpositions`<- lapply(new_positive$`BT_processed_positive$totalpositions`, gsub, pattern=c(" "), replacement='')

write(unlist(new_positive$`BT_processed_positive$totalpositions`), file = "positive_positions.fasta", sep = ">")

#position mentioned in positive_positions.fasta for sequence chr1_5583493_5588992 is not valid. Letter 'C' was not found at position 5161

# checking where the C's are
chr1_5583493_5588992 <- gregexpr(pattern ='C', BT_processed_positive$seq[[1]])
chr1_5583493_5588992_t <- tibble(unlist(chr1_5583493_5588992))
#there is a c at 5162!


##### retry with only negative strands #####

write.fasta( as.list(BT_processed_negative$seq), BT_processed_negative$names, file.out = "negative_seq_only.fasta", as.string = T)

modified_BT_processed_negative <- tibble(BT_processed_negative$totalpositions) %>% 
  add_column(BT_processed_negative$names)

# double all rows 
new_negative <- modified_BT_processed_negative[rep(1:nrow(modified_BT_processed_negative),1,each=2),]

# replace all duplicates with appropriate sequence names
new_negative$`BT_processed_negative$totalpositions`[c(seq(1, dim(new_negative)[1], by=2))] <- paste(">", new_negative$`BT_processed_negative$names`[c(seq(2, dim(new)[1], by=2))])
new_negative$`BT_processed_negative$totalpositions`<- lapply(new_negative$`BT_processed_negative$totalpositions`, gsub, pattern=c(" "), replacement='')

write(unlist(new_negative$`BT_processed_negative$totalpositions`), file = "negative_positions.fasta", sep = ">")

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

##### DNASHAPER #####

pred_methy <- getShape("positive_seq_only.fasta", methylate = TRUE, methylatedPosFile = "positive_positions.fasta")
plotShape(pred_methy$MGW)

pred_methy_neg <- getShape("negative_seq_only.fasta", methylate = TRUE, methylatedPosFile = "negative_positions.fasta")
plotShape(pred_methy_neg$MGW)

# MGW can be replaced with HelT, Rise, Roll, Shift, Slide, Tilt, Buckle, Opening, ProT, Shear, Stagger, Stretch or EP
plotShape(pred_methy$HelT)
plotShape(pred_methy$ProT)
plotShape(pred_methy$Roll)

# plotShape(pred_methy$Tilt) ### Error in colMeans(shapeMatrix, na.rm = TRUE) : 
# 'x' must be an array of at least two dimensions
# In addition: There were 15 warnings (use warnings() to see them)

heatShape(pred_methy$ProT, 10)

### will need to change the data shape before doing this but I think it will be a good diagram
# *Note that the input data should only contain one sequence.
# fn2 <- system.file("extdata", "SingleSeqsample.fa", package = "DNAshapeR")
# pred2 <- getShape(fn2)
# trackShape(fn2, pred2) # Only for single sequence files




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