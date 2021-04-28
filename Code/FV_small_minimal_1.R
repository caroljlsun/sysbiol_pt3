# TO DO BEFORE START

# save feature vectors from every piece of script I have written. I will probably make a few more than I already have for all the B, T cells, and the multiple features that DNAshaper has

#probably need to do a rolling window for GC content, cuz it makes no sense to do it 1-mer
#not sure how big the window should be, could go for 30 cuz it sounds right
#and make changes accordingly

#after I've saved a million smaller feature vectors

#import them here, and we will construct the most unholy behemoth of feature vectors

# load workspace
setwd("/home/cjls4/ML/")
load("attempt_1")
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
library(stringr)

library(learnMotifs)
library(mltools)

library(caret)
library(keras)

library(reticulate)
library(tensorflow)
library(useful)

#Global options

setwd("/home/cjls4/feature_vectors/")

#uhhh

install.packages("kerasR")

devtools::install_github("rstudio/reticulate")
library(reticulate)
#use_python("/home/cjls4/miniconda3/envs/r-reticulate/bin/python3")
# py_config()
# 
# library(reticulate)
use_python("/usr/bin/python3")
# py_config()
# 
# use_condaenv("tf")
# tensorflow::install_tensorflow("auto")
# reticulate::py_discover_config()

tf$constant("Hellow Tensorflow")

#I will probably need to drag in the tibbles instead of just the fv's lol

#Let's do a toy example and make a FV with both the expression data and the G4 quadruplexes.
#Somehow need to match the sequences/keys before making the actual FV
# I think it should look something like this
# A T C G BF BM TF TM G4
# 1 0 0 0 2  0  0  1  0
# 0 0 0 1 1  1  1  1  0.1
# 0 0 1 0 0  0  1  1  0.2

#with a new page? whatever, per promoter
# so I should end up with a dataset of fvs
# which have dimensions of 58879x2749x9
# or more aesthetically, 9*2749*58879
# and each time I add more to this FV, it should only expand the columns
# so, n*2749*58879

#remember to exclude chromsomes 1 and 2. 1 for test, 2 for validation or something

#import data

load("minus_k_G4_final")

load("expression_tibble.RData")

#uhh <- minus_k_fv_tibble[minus_k_fv_tibble$chromosome %in% intersect_tibble$names]

#G4_exp <- inner_join(intersect_tibble, minus_k_fv_tibble, by = c("names" = "chromosome"))

G4_exp_join <- inner_join(intersect_tibble, minus_k_fv_final, by = c("seq" = "seq"))
G4_exp_join <- G4_exp_join %>% 
  dplyr::select(-chromosome)

length(unique(G4_exp_join$seq)) #5748

hold <- tibble(intersect(intersect_tibble$seq, minus_k_fv_final$seq))


##### Learnmotifs has the best method for one hot encoding #####
sampled <- G4_exp_join[1,]

s1 <- learnMotifs::one_hot(sampled$seq, zeros_len = 5)

s2 <- as.data.table( t(s1[,1:4, 5:2754,]))

#make a column of ATGC

s3 <- unlist(strsplit(sampled$seq, split = ""))
s3 <- transpose(as.data.table(s3))

#make a collection of info, no padding. We can add padding later if need be

sampled_fv <- tibble(sampled) %>% 
  dplyr::select(-names, -seq, -bcellf, -bcellm, -tcellf, -tcellm) %>%
  t()%>% 
  as.data.table() %>% 
  #rename(G4 = V1) %>% 
  add_column(s2, .before = 2) %>% 
  #add_column(s3, .before = 1) %>% 
  add_column(rep(sampled$bcellf, times = 2750)) %>% 
  add_column(rep(sampled$bcellm, times = 2750)) %>% 
  add_column(rep(sampled$tcellf, times = 2750)) %>% 
  add_column(rep(sampled$tcellm, times = 2750)) %>% 
  dplyr::rename(seq = V1.1,
         bcellf = `rep(sampled$bcellf, times = 2750)`,
         bcellm = `rep(sampled$bcellm, times = 2750)`,
         tcellf = `rep(sampled$tcellf, times = 2750)`,
         tcellm = `rep(sampled$tcellm, times = 2750)`)

#Do the above again, but now on all of the data...

#let's try it with a 10 set of data, to practice the iterative steps etc.
# from last time, the best way to do things was use a dataframe, and set everything to numeric if possible

biggersample <- G4_exp_join[1:10,]

#something is wrong with the expression data wtf
#actually it's fine dw

bigger_fv <- data.table()

#as.data.table( t(s1[,1:4, 5:2754,])
               
custom_one_hot <- function(i) {
bigger_fv[i] <- as.data.table(t(learnMotifs::one_hot(biggersample$seq[i], zeros_len = 5)[,1:4, 5:2754,]))  
}

b1 <- sapply(1:10, custom_one_hot)
b1.1 <- as.data.table(b1)


# really think about what shape this will be
# 
# w1 <- simplify2array(biggersample)
# w2 <- base::split(biggersample, biggersample$seq)

biggersample_no_exp <- biggersample %>% 
  dplyr::select(-bcellf, -bcellm, -tcellf, -tcellm, -seq, -names) %>% 
  transpose()

#w3 <- abind(base::split(biggersample_no_exp, biggersample_no_exp$seq), along = 3)

biggersample_exp <- tibble(rep(biggersample$bcellf, times = 2750))

w4 <- array(NA, c(2750,9,10))

w4[1:2750,1,1:10] <- unlist(biggersample_no_exp)

what <- rep(biggersample$bcellf, each = 2750)

w4[1:2750,2,1:10] <- what

w5 <- w4
# this doesn't work, it adds rowwise instead of colwise
# more_expression <- unlist(c(rep(biggersample$bcellm, each = 2750),
#                             rep(biggersample$tcellf, each = 2750),
#                             rep(biggersample$tcellm, each = 2750)))
# 
# w5[1:2750,3:5, 1:10] <- more_expression

w4[1:2750,3,1:10] <- rep(biggersample$bcellm, each = 2750)
w4[1:2750,4,1:10] <- rep(biggersample$tcellf, each = 2750)
w4[1:2750,5,1:10] <- rep(biggersample$tcellm, each = 2750)

w6 <- w4

w6[1:2750,6,1:10] <- unlist(b1.1[1,])
w6[1:2750,7,1:10] <- unlist(b1.1[2,])
w6[1:2750,8,1:10] <- unlist(b1.1[3,])
w6[1:2750,9,1:10] <- unlist(b1.1[4,])


# What I did previously, adapted


set.seed(12345)

sample_train <- w6[,,4:10]
sample_test <- w6[,,1:3]

x_train <-sample_train[,2:9,1:7]
y_train <-sample_train[1:2750,1,1:7]

x_uhh_train <- as.matrix(x_train)
y_uhh_train <- as.matrix(y_train)

x_test <- sample_test[,2:9,]
y_test <- sample_test[,1,]

#dimensions declarations

dim(x_train)#2750,8,7
dim(y_train)
#keras isn't very happy with the dimensions I have right now
# I suspect I need to change it to matrix slice, row, column

#and it is currently row, col, matrix slice. 

please_x_train <- aperm(x_train, c(3,1,2))
dim(please_x_train)

#normalise data?

norm_please_x_train <-  normalize(please_x_train, axis = -1, order = 2)

#export some of this data, using in attempt_keras_2

save(x_train, file = "x_train.RData")
save(x_test, file = "x_test.RData")

save(y_train, file = "y_train.RData")
save(y_test, file = "y_test.RData")


#keras troubleshooting


setwd("/home/cjls4/ML/")
#save.image("attempt_1")

#Ok make a full size dataset



