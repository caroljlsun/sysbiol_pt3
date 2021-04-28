

# load workspace
setwd("/home/cjls4/ML/FV_bigger_minimal_1/")
load("FV_bigger_minimal_1.RData")

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
library(abind)

library(learnMotifs)
library(mltools)

library(keras)

library(reticulate)
library(tensorflow)
library(useful)

#Global options

setwd("/home/cjls4/feature_vectors/")

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

G4_exp_antijoin <- anti_join(intersect_tibble, minus_k_fv_final, by = c("seq" = "seq"))


length(unique(G4_exp_join$seq)) #5748

hold <- tibble(intersect(intersect_tibble$seq, minus_k_fv_final$seq))


#Let's get an array for all the G4 positive G4s

full_minimal <- G4_exp_join

full_minimal_fv <- data.table()

fm_one_hot <- function(i) {
  full_minimal_fv[i] <- as.data.table(t(learnMotifs::one_hot(full_minimal$seq[i], zeros_len = 5)[,1:4, 5:2754,]))  
}

fm1 <- sapply(1:dim(full_minimal)[1], fm_one_hot)

fm1.1 <- as.data.table(fm1)

#p1 = positive G4s 
#will have 8 columns, 4 for the expression, 4 for the hot encoding
p1 <- array(NA, c(2750,8,5748))

p1[1:2750,1,1:5748] <- unlist(fm1.1[1,]) #A
p1[1:2750,2,1:5748] <- unlist(fm1.1[2,]) #C
p1[1:2750,3,1:5748] <- unlist(fm1.1[3,]) #G
p1[1:2750,4,1:5748] <- unlist(fm1.1[4,]) #T

p1[1:2750,5,1:5748] <- rep(full_minimal$bcellf, each = 2750)

p1[1:2750,6,1:5748] <- rep(full_minimal$bcellm, each = 2750)

p1[1:2750,7,1:5748] <- rep(full_minimal$tcellf, each = 2750)

p1[1:2750,8,1:5748] <- rep(full_minimal$tcellm, each = 2750)


# make an array with all promoters that DON'T have G4s

negative_minimal <- G4_exp_antijoin

negative_minimal_fv <- data.table()

nm_one_hot <- function(i) {
  negative_minimal_fv[i] <- as.data.table(t(learnMotifs::one_hot(negative_minimal$seq[i], zeros_len = 5)[,1:4, 5:2754,]))  
}

nm1 <- sapply(1:dim(negative_minimal)[1], nm_one_hot)

nm1.1 <- as.data.table(nm1)

#n1 = negative G4s 
#will have 8 columns, 4 for the expression, 4 for the hot encoding#

#btw check this negative set for duplicates
#length(unique(G4_exp_antijoin$seq))

n1 <- array(NA, c(2750,8,17516))

n1[1:2750,1,1:17516] <- unlist(nm1.1[1,]) #A
n1[1:2750,2,1:17516] <- unlist(nm1.1[2,]) #C
n1[1:2750,3,1:17516] <- unlist(nm1.1[3,]) #G
n1[1:2750,4,1:17516] <- unlist(nm1.1[4,]) #T


browseURL('https://www.youtube.com/watch?v=QH2-TGUlwu4')

n1[1:2750,5,1:17516] <- rep(negative_minimal$bcellf, each = 2750)

n1[1:2750,6,1:17516] <- rep(negative_minimal$bcellm, each = 2750)

n1[1:2750,7,1:17516] <- rep(negative_minimal$tcellf, each = 2750)

n1[1:2750,8,1:17516] <- rep(negative_minimal$tcellm, each = 2750)


# add the negative G4 array to the positive

pn1 <- abind(p1, n1)
#god it's massive lol

#exclude chromosome 1 from the training

p_chr1s <- str_detect(full_minimal$names, "chr1_")

#I think I need the opposite to subset the fv with
p_chr1s_opp <- str_detect(full_minimal$names, "chr1_", negate = TRUE)

table(p_chr1s)
table(p_chr1s_opp)


n_chr1s <- str_detect(negative_minimal$names, "chr1_")
n_chr1s_opp <- str_detect(negative_minimal$names, "chr1_", negate = TRUE)

table(n_chr1s)
table(n_chr1s_opp)

# time to pull these promoters out

pn_chr1s_opp <- c(p_chr1s_opp, n_chr1s_opp)

x_train <- pn1[1:2750, 1:8, pn_chr1s_opp]
dim(x_train)

y_train <-  c(rep(1, times = 5436), rep(0, times = 16477))


pn_chr1s <- c(p_chr1s, n_chr1s)

x_test <- pn1[1:2750, 1:8, pn_chr1s]
dim(x_test)

y_test <- c(rep(1, times = 312 ), rep(0, times = 1039))


# save into x y train test for keras
# will be use in attempt_keras_3

setwd("/home/cjls4/ML/FV_bigger_minimal_1/")
save(x_train, file = "x_train.RData")

save(y_train, file = "y_train.RData")

save(x_test, file = "x_test.RData")

save(y_test, file = "y_test.RData")

#save everything

save.image(file = "FV_bigger_minimal_1.RData")

#I will probably take a few more chromosomes out, hmmm

p_chr2s <- str_detect(full_minimal$names, "chr2_")

table(p_chr2s)
# maybe later, I just want to see if it works



#do a quick tally of number of promoters per chromosome

c1s <- str_detect(full_minimal$names, "chr1_")
c19s <- str_detect(full_minimal$names, "chr19_")
cXs <- str_detect(full_minimal$names, "chrX_")
cYs <- str_detect(full_minimal$names, "chrY_")

cs <- vector()
tcs <- vector()
for (i in 1:19) {
  cs[i] <- list(str_detect(full_minimal$names, paste0("chr", i, "_")))
  tcs[i] <- table(cs[i])[[2]] 
}

table(cs[1])


table(cs[1])[[2]]


tcs[20] <- table(cXs)[[2]]

tcs <- as.data.table(tcs)

chr_counts <- ggplot(data = tcs, mapping = aes(x = 1:20, y=tcs ))+
  geom_col()
chr_counts


#now for the antijoin

ncs <- vector()
ntcs <- vector()
for (i in 1:19) {
  ncs[i] <- list(str_detect(negative_minimal$names, paste0("chr", i, "_")))
  ntcs[i] <- table(ncs[i])[[2]] 
}

ncXs <- str_detect(negative_minimal$names, "chrX_")
ncYs <- str_detect(negative_minimal$names, "chrY_")

table(ncYs)[[2]]



ntcs[20] <- table(ncXs)[[2]]
ntcs[21] <- table(ncYs)[[2]]


ntcs <- as.data.table(ntcs)

n_chr_counts <- ggplot(data = ntcs, mapping = aes(x = 1:21, y=ntcs ))+
  geom_col()
n_chr_counts

#both?

add_tcs <- c(unlist(tcs), 0)
add_tcs <- append(add_tcs, as.integer(0))

uhh <- cbind(add_tcs, ntcs)
then <- rowSums(uhh)

then <- as.data.frame(then)

t_chr_counts <- ggplot(data = then, mapping = aes(x = 1:21, y = then))+
  geom_col()
t_chr_counts
