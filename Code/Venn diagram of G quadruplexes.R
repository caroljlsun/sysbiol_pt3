load("Venn.RData")

#libraries

library(dplyr)
library(useful)
library(ggplot2)
library(ggvenn)
library(ggVennDiagram)
library(viridis)
library(tidyverse)
library(plyr)
library(stringr)
library(seqinr)
library(Biostrings)
library(data.table)
#Global options

setwd(dir = "/home/cjls4/feature_vectors/")
options(scipen=999)
options(expressions = 500000)

#Import data

load("minus_k_G4_final")
load("plus_k_G4_final")

load("minus_PDS_G4_final")
load("plus_PDS_G4_final")

##### how much do the minus and plus k's overlap? #####
#may need to use reverse complement before I can do that

#h1 <- reverseComplement(DNAString(x=minus_k_fv_final$seq[1], start = 1))
#h2 <- data.frame(as.character(h1))

reverse_plus_k_fv_final <- data.frame(matrix(data = NA, nrow = dim(plus_k_fv_final)[1], ncol = 1))

for (i in 1:dim(plus_k_fv_final)[1]){
  reverse_plus_k_fv_final[i,] <- as.character(reverseComplement(DNAString(x = plus_k_fv_final$seq[i])))
}

names(reverse_plus_k_fv_final) <- "seq"

intersecting_f_r_k <- inner_join(minus_k_fv_final, reverse_plus_k_fv_final, by = c("seq" = "seq"))

#I guess they don't have many, literally one?


##### how much do the minus and plus PDS's overlap?#####

reverse_plus_PDS_fv_final <- data.frame(matrix(data = NA, nrow = dim(plus_PDS_fv_final)[1], ncol = 1))

for (i in 1:dim(plus_PDS_fv_final)[1]) {
  reverse_plus_PDS_fv_final[i,] <-as.character(reverseComplement(DNAString(x = plus_PDS_fv_final$seq[i])))
}

names(reverse_plus_PDS_fv_final) <- "seq"

intersecting_f_r_PDS <- inner_join(minus_PDS_fv_final, reverse_plus_PDS_fv_final, by = c("seq" = "seq"))

#why is it the same one as the K ones

##### how much do the minus k's and PDS's overlap? #####

intersecting_m_k_pds <- inner_join(minus_k_fv_final, minus_PDS_fv_final, by = c("seq" = "seq"))

#double check the intersect with another method

int_1 <- intersect(minus_k_fv_final$seq, minus_PDS_fv_final$seq)

int_2 <- data.frame(int_1)

#and then find out how many belong to K and PDS, respectively


intersect_1 <- intersect(minus_k_fv_final$seq, intersecting_m_k_pds$seq)
dim(intersecting_m_k_pds)
dim(minus_k_fv_final)

# 355 promoters are from PDS, 10062 are shared between minus K and PDS

intersect_2 <- intersect(minus_PDS_fv_final$seq, intersecting_m_k_pds$seq)
dim(minus_PDS_fv_final)


##### Save the minus intersecting G4s #####

minus_intersect_G4s <- intersecting_m_k_pds[,1:2752]

save(minus_intersect_G4s, file = "minus_intersect_G4s.RData")


d1 <- t(minus_intersect_G4s) %>% 
  as.data.table()

d2 <- minus_intersect_G4s[,3:2752] 

d3 <- transpose(d2)


d3$position <- seq(from = -2500, to = 249, by = 1)

d3$average <- apply(d3[1:10062], 1, mean)

library(ggplot2)

p1 <- ggplot(data = d3, mapping = aes(y = average, x = position)) +
  geom_line()


p1

##### how much do the plus k's and PDS's overlap? #####

intersecting_p_k_pds <- inner_join(plus_k_fv_final, plus_PDS_fv_final, by = c("seq" = "seq"))

# how many belong to K and PDS?

intersect_3 <- intersect(plus_k_fv_final$seq, intersecting_p_k_pds$seq)
length(intersect_3) #9978
dim(intersecting_p_k_pds) #9978
dim(plus_k_fv_final)#10336

intersect_4 <- intersect(plus_PDS_fv_final$seq, intersecting_p_k_pds$seq)
length(intersect_4) #9978
dim(plus_PDS_fv_final)#16161

##### Save the plus intersecting G4s #####

plus_intersect_G4s <- intersecting_p_k_pds[, 1:2752]

save(plus_intersect_G4s, file = "plus_intersect_G4s.RData")

save.image(file = "Venn.RData")
