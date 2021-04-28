library(GenomicRanges)
library(genomation)
library(readr)
library(BiSeq)

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

sample_id_1 <- "sample_1"
treatment_1 <- "treatment_B_F"

#reassign names

rename_methyl <- function(object)
{
  propernames <- list("chr", "start", "end", "strand", "coverage", "numCs", "numTs")
  
  for(i in 1:length(propernames)){
    object@names[i] <- propernames[[i]]
  }
  return(object)
}



#just to see if it is working
why <- readBismarkProcessed("B6_F_B_BS.10.methy.combined_strand.mincov_10.tsv",
                            sample.id =  sample_id_1,
                            treatment =  treatment_1)

why2 <- rename_methyl(why)

getMethylationStats(why, plot = T)
getCoverageStats(why, plot = T)


# ok let's get several .tsv's converted

sample_ids_2 <- list("sample_F_B_10", "sample_F_B_12", "sample_F_T_13", "sample_F_T_15")
treatments_2 <- c(1,1,0,0)

small_group_of_tsv <- list("B6_F_B_BS.10.methy.combined_strand.mincov_10.tsv",
                           "B6_F_B_BS.12.methy.combined_strand.mincov_10.tsv",
                           "B6_F_T_BS.13.methy.combined_strand.mincov_10.tsv",
                           "B6_F_T_BS.15.methy.combined_strand.mincov_10.tsv")

small_group_methyl <- read_bismark_processed(small_group_of_tsv,
                                           sample.id = sample_ids_2,
                                           treatment = treatments_2)

# renaming the header in all of the samples

for (i in 1:length(small_group_methyl)){
  small_group_methyl[[i]] <- rename_methyl(small_group_methyl[[i]])
}

small_group_methyl <- as(small_group_methyl, "methylRawList")

small_group_unite <- unite(small_group_methyl, save.db = T)

#small_group_methyl_DB <- as(small_group_methyl, "methylRawListDB")

getMethylationStats(small_group_methyl[[4]], plot = T)

getCorrelation(small_group_unite, plot = T)

clusterSamples(small_group_unite, dist = "correlation", method = "ward", plot = T)

PCASamples(small_group_unite)

# first reduce/compile into B and T cells



# then reduce/compile into FB, FT, MB, MT



#correlate with genomic positions (specifically the promoter regions)

# need to get a .bed file?

annotateWithFeature()

annotateWithGeneParts()
# HMMMMMMM