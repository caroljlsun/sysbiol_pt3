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
library(plyr)

library(learnMotifs)
library(mltools)


#Global options

setwd("/home/cjls4/feature_vectors/")

#load("Minus_G4_positive.RData")

# I think it should look something like this
# A T C G BF BM TF TM H3Kac H3Kac H3Kme H3kme RL M1 M2 M3 M4 Methyl GC  
# 1 0 0 0 2  0  0  1  0     2     2     3     1 2.3  2 0  0  0      0
# 0 0 0 1 1  1  1  1  0.1   3     1     0     1  2   2 2  3  0      1
# 0 0 1 0 0  0  1  1  0.2   0     0     0     0  1   2 2  2  3      1


#remember to exclude chromsomes 1 and 2. 1 for test, 2 for validation or something

#import data needed for the minuses

load("minus_intersect_G4s.RData")
load("expression_tibble.RData")
load("H3K27ac_B_final")
load("H3K27ac_T_final")
load("H3K27me3_B_final")
load("H3K27me3_T_final")
load("loops_fv_final")

load("vector_final_pos_MGW")
load("vector_final_pos_ProT")
load("vector_final_pos_Roll")
load("vector_final_pos_HelT")

load("methyl_fv_final")
load("minus_intersect_gc")
load("plus_intersect_gc")

##### Main G4 presence data #####

G4_exp_join_min <- inner_join(intersect_tibble, minus_intersect_G4s, by = c("seq" = "seq"))

G4_exp_join_min <- G4_exp_join_min %>% 
  dplyr::select(-chromosome.x)

#Let's get an array for all the G4 positive G4s

full_minimal <- G4_exp_join_min

full_minimal_fv <- data.table()

##### one hot encoding + expression #####

fm_one_hot <- function(i) {
  full_minimal_fv[i] <- as.data.table(t(learnMotifs::one_hot(full_minimal$seq[i], zeros_len = 5)[,1:4, 5:2754,]))  
}

fm1 <- sapply(1:dim(full_minimal)[1], fm_one_hot)

fm1.1 <- as.data.table(fm1)

##### Histone marks #####

check1 <- inner_join(G4_exp_join_min, H3K27ac_B_final, by = c("seq" = "seq"))

#h3k27ac_b

a1 <- base::merge(G4_exp_join_min, H3K27ac_B_final, "seq", all.x = T)
a2 <- a1 %>% 
  dplyr::select(V1:V2750)
a2[is.na(a2)] <- 0

a3 <- array(NA, c(5555,2750,1))
a3[,,1] <- unlist(a2)

a4 <- aperm(a3, c(2,3,1))

# h3k27ac_t

b1 <- base::merge(G4_exp_join_min, H3K27ac_T_final, "seq", all.x = T)
b2 <- b1 %>% 
  dplyr::select(V1:V2750)
b2[is.na(b2)] <- 0

b3 <- array(NA, c(5555,2750,1))
b3[,,1] <- unlist(b2)

b4 <- aperm(b3, c(2,3,1))

# h3k27me3_b

c1 <- base::merge(G4_exp_join_min, H3K27me3_B_final, "seq", all.x = T)
c2 <- c1 %>% 
  dplyr::select(V1:V2750)
c2[is.na(c2)] <- 0

c3 <- array(NA, c(5555,2750,1))
c3[,,1] <- unlist(c2)

c4 <- aperm(c3, c(2,3,1))

# h3k27me3_t

d1 <- base::merge(G4_exp_join_min, H3K27me3_T_final , "seq", all.x = T)
d2 <- d1 %>% 
  dplyr::select(V1:V2750)
d2[is.na(d2)] <- 0

d3 <- array(NA, c(5555,2750,1))
d3[,,1] <- unlist(d2)

d4 <- aperm(d3, c(2,3,1))

##### R loops #####

r1 <- base::merge(G4_exp_join_min, loops_fv_final, "seq", all.x = T)
r2 <- r1 %>% 
  dplyr::select(V1:V2750)
r2[is.na(r2)] <- 0

r3 <- array(NA, c(5555,2750,1))
r3[,,1] <- unlist(r2)

r4 <- aperm(r3, c(2,3,1))

##### DNAshapeR #####

check3 <- inner_join(G4_exp_join_min, vector_final_pos_MGW, by = c("seq" = "seq"))

blank_columns <- as.data.frame(matrix(data = NA, nrow = 5555, ncol = 2))

blank_column_1 <- as.data.frame(matrix(data = NA, nrow = 5555, ncol = 1))

#MGW

MGW1m <- base::merge(G4_exp_join_min, vector_final_pos_MGW , "seq", all.x = T)
MGW2m <- MGW1m %>% 
  dplyr::select(V1:V2746) %>%
  add_column(blank_columns, .before = 1) %>% 
  add_column(blank_columns)


MGW2m[is.na(MGW2m)] <- 0

dim(MGW2m) # 5555 2750

MGW2m_1 <- array(NA, c(5555,2750,1))

MGW2m_1[,,1] <- unlist(MGW2m)

MGW2m_2 <- aperm(MGW2m_1, c(2,3,1))

#ProT

ProT1m <- base::merge(G4_exp_join_min, vector_final_pos_ProT , "seq", all.x = T)

ProT2m <- ProT1m %>% 
  dplyr::select(V1:V2746) %>%
  add_column(blank_columns, .before = 1) %>% 
  add_column(blank_columns)

ProT2m[is.na(ProT2m)] <- 0

dim(ProT2m) # 5555 2750

ProT2m_1 <- array(NA, c(5555,2750,1))

ProT2m_1[,,1] <- unlist(ProT2m)

ProT2m_2 <- aperm(ProT2m_1, c(2,3,1))

#Roll
#Since Roll is one of the base pair step features, it's size is kind of off. Have decided to use zero padding
# one blank column at the start, two at the end.

Roll1m <- base::merge(G4_exp_join_min, vector_final_pos_Roll , "seq", all.x = T)

Roll2m <- Roll1m %>% 
  dplyr::select(V1:V2747) %>%
  add_column(blank_column_1, .before = 1) %>% 
  add_column(blank_columns)

Roll2m[is.na(Roll2m)] <- 0

dim(Roll2m) # 5555 2750

Roll2m_1 <- array(NA, c(5555,2750,1))

Roll2m_1[,,1] <- unlist(Roll2m)

Roll2m_2 <- aperm(Roll2m_1, c(2,3,1))

#HelT
# same as above

HelT1m <- base::merge(G4_exp_join_min, vector_final_pos_HelT , "seq", all.x = T)

HelT2m <- HelT1m %>% 
  dplyr::select(V1:V2747) %>%
  add_column(blank_column_1, .before = 1) %>% 
  add_column(blank_columns)

HelT2m[is.na(HelT2m)] <- 0

dim(HelT2m) # 5555 2750

HelT2m_1 <- array(NA, c(5555,2750,1))

HelT2m_1[,,1] <- unlist(HelT2m)

HelT2m_2 <- aperm(HelT2m_1, c(2,3,1))

##### Methylated positions #####

check3 <- inner_join(G4_exp_join_min, methyl_fv_final, by = c("seq" = "seq"))

m1 <- base::merge(G4_exp_join_min, methyl_fv_final, "seq", all.x = T)
m2 <- m1 %>% 
  dplyr::select(V1:V2750)
m2[is.na(m2)] <- 0

dim(m2) # 5555 2750

m3 <- array(NA, c(5555,2750,1))
m3[,,1] <- unlist(m2)

m4 <- aperm(m3, c(2,3,1))

##### GC content #####

check4 <- inner_join(G4_exp_join_min, minus_intersect_gc, by = c("seq" = "seq"))

gc1 <- base::merge(G4_exp_join_min, minus_intersect_gc, "seq", all.x = T)
gc2 <- gc1 %>% 
  dplyr::select(V1:V110)

gc3 <- gc2 %>% 
  rep(each = 25)

gc4 <- as.data.table(gc3)

#just checking where the NAs come from. It's from all the NNNNATGCGCT... sequences. Will be assigned as 0 in the final
# assembled feature vector
#hmmmmm <- is.na(gc4)
#uhhhhh <- which(is.na(gc2))

gc4.1 <- ldply(gc3)

dim(gc4) # 5555 2750
# row, column, slice
gc5 <- array(NA, c(5555,2750,1))
gc5[,,1] <- unlist(gc4)

gc6 <- aperm(gc5, c(2,3,1))


##### The actual FV ##### 
#min_p1 = positive G4s, minus strand sequences
#will have 19 columns
min_p1 <- array(NA, c(2750,19,5555))

min_p1[1:2750,1,1:5555] <- unlist(fm1.1[1,]) #A
min_p1[1:2750,2,1:5555] <- unlist(fm1.1[2,]) #C
min_p1[1:2750,3,1:5555] <- unlist(fm1.1[3,]) #G
min_p1[1:2750,4,1:5555] <- unlist(fm1.1[4,]) #T

min_p1[1:2750,5,1:5555] <- rep(full_minimal$bcellf, each = 2750) #bcellf expression
min_p1[1:2750,6,1:5555] <- rep(full_minimal$bcellm, each = 2750) #bcellm expression
min_p1[1:2750,7,1:5555] <- rep(full_minimal$tcellf, each = 2750) #tcellf expression
min_p1[1:2750,8,1:5555] <- rep(full_minimal$tcellm, each = 2750) #tcellm expression

min_p1[1:2750,9,1:5555] <- a4
min_p1[1:2750,10,1:5555] <- b4
min_p1[1:2750,11,1:5555] <- c4
min_p1[1:2750,12,1:5555] <- d4

min_p1[1:2750,13,1:5555] <- r4

min_p1[1:2750,14,1:5555] <- MGW2m_2
min_p1[1:2750,15,1:5555] <- ProT2m_2
min_p1[1:2750,16,1:5555] <- Roll2m_2
min_p1[1:2750,17,1:5555] <- HelT2m_2

min_p1[1:2750,18,1:5555] <- m4

min_p1[1:2750,19,1:5555] <- gc6


min_p1[1:100,,1]

# exclude chromosome 1 from the training

p_chr1s <- str_detect(full_minimal$names, "chr1_")

p_chr1s_opp <- str_detect(full_minimal$names, "chr1_", negate = TRUE)

table(p_chr1s)
table(p_chr1s_opp)

# no chr1s
x_train_min_p <- min_p1[1:2750, 1:19, p_chr1s_opp]

y_train_min_p <-  c(rep(1, times = 5257))


# only chr1s
x_test_min_p <- min_p1[1:2750, 1:19, p_chr1s]

y_test_min_p <- c(rep(1, times = 298 ))


save(x_train_min_p, file = "x_train_min_p")
save(y_train_min_p, file = "y_train_min_p")
save(x_test_min_p, file = "x_test_min_p")
save(y_test_min_p, file = "y_test_min_p")


save(min_p1, file = "MINUS_STRAND_G4_POSITIVE.RData")
save.image("Minus_G4_positive.RData")



##### get sequences separately #####

MINUS_STRAND_G4_POSITIVE_SEQ <- G4_exp_join_min$seq %>% 
  as.data.frame()

MINUS_STRAND_G4_POSITIVE_SEQ$row_names <- row.names(MINUS_STRAND_G4_POSITIVE_SEQ)

MINUS_STRAND_G4_POSITIVE_SEQ <- tibble(MINUS_STRAND_G4_POSITIVE_SEQ) %>% 
  dplyr::rename(.,seq = `.`)

MINUS_STRAND_G4_POSITIVE_SEQ$seq <- as.character(MINUS_STRAND_G4_POSITIVE_SEQ$seq)

MINUS_STRAND_G4_POSITIVE_SEQ <- as.data.frame(MINUS_STRAND_G4_POSITIVE_SEQ)

MINUS_STRAND_G4_POSITIVE_SEQ$G4 <- 1

MINUS_STRAND_G4_POSITIVE_SEQ$strand <- "-"

save(MINUS_STRAND_G4_POSITIVE_SEQ, file = "MINUS_STRAND_G4_POSITIVE_SEQ.RData")

##### get the gene names separately #####

MINUS_STRAND_G4_POSITIVE_NAMES <- G4_exp_join_min$gene_names %>% 
  as.data.table

names(MINUS_STRAND_G4_POSITIVE_NAMES) <- "gene_names"

save(MINUS_STRAND_G4_POSITIVE_NAMES, file = "MINUS_STRAND_G4_POSITIVE_NAMES.RData")


min_g4_positive_names_train <- MINUS_STRAND_G4_POSITIVE_NAMES[p_chr1s_opp]
save(min_g4_positive_names_train, file = "min_g4_positive_names_train.RData")

min_g4_positive_names_test <- MINUS_STRAND_G4_POSITIVE_NAMES[p_chr1s]
save(min_g4_positive_names_test, file = "min_g4_positive_names_test.RData")


