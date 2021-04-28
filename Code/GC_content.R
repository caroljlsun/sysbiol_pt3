# have a look for packages which calc gc for you

#look into using perl
# in gsub

load("GC.RData")

library(seqinr)
library(useful)
library(dplyr)

setwd("/home/cjls4/feature_vectors/")

load("minus_intersect_G4s.RData")
load("expression_tibble.RData")

# minus strand, presence of G4

window_size <- 25

chunks <- seq(from = 1, to = 2750, by = window_size)

g1 <- matrix(data = 0, nrow = dim(minus_intersect_G4s)[1], ncol = 110)


for (j in 1:dim(minus_intersect_G4s)[1]) {
  for (i in 1:110) {
    g1[j,i] <- GC(s2c(minus_intersect_G4s$seq[j])[(chunks[i]):(chunks[i]+window_size-1)])
  }
  g1
}

# replace NAs

g1[is.na(g1)] <- 0


# assign to the right chromosome/sequences

minus_intersect_gc <- tibble(as.data.table(g1)) %>% 
  add_column(minus_intersect_G4s$seq, .before = 1) %>% 
  add_column(minus_intersect_G4s$chromosome.x, .before = 1) %>% 
  dplyr::rename(chromosome = `minus_intersect_G4s$chromosome.x`) %>% 
  dplyr::rename(seq = `minus_intersect_G4s$seq`)

#find and eliminate the NAs

hm1 <- which(is.na(g1))
# the mystery has been solved! the sequences sometimes are NNNNNTGCCGTAA... meaning there is no GC score
# therefore just replace the NAs with 0


#save the gc fv

save(minus_intersect_gc, file = "minus_intersect_gc")

#minus strand, absence of G4

G4_exp_join_absence <- anti_join(intersect_tibble, minus_intersect_G4s, by = c("seq" = "seq"))

g1_1 <- matrix(data = NA, nrow = dim(G4_exp_join_absence)[1], ncol = 110)


for (j in 1:dim(G4_exp_join_absence)[1]) {
  for (i in 1:110) {
    g1_1[j,i] <- GC(s2c(G4_exp_join_absence$seq[j])[(chunks[i]):(chunks[i]+window_size-1)])
  }
  g1_1
}

# replace NAs

g1_1[is.na(g1_1)] <- 0


# assign to the right chromosome/sequences

minus_intersect_absence_gc <- tibble(as.data.table(g1_1)) %>% 
  add_column(G4_exp_join_absence$seq, .before = 1) %>% 
  add_column(G4_exp_join_absence$names, .before = 1) %>% 
  dplyr::rename(chromosome = `G4_exp_join_absence$names`) %>% 
  dplyr::rename(seq = `G4_exp_join_absence$seq`)

#save

save(minus_intersect_absence_gc, file = "minus_intersect_absence_gc.RData")


# plus strand, presence of G4s

load("plus_intersect_G4s.RData")

g2 <- matrix(data = NA, nrow = dim(plus_intersect_G4s)[1], ncol = 110)

for (j in 1:dim(plus_intersect_G4s)[1]) {
  for (i in 1:110) {
    g2[j,i] <- GC(s2c(plus_intersect_G4s$seq[j])[(chunks[i]):(chunks[i]+window_size-1)])
  }
  g2
}

# replace NAs

g2[is.na(g2)] <- 0


# assign to the right chromosome/sequences

plus_intersect_gc <- tibble(as.data.table(g2)) %>% 
  add_column(plus_intersect_G4s$seq, .before = 1) %>% 
  add_column(plus_intersect_G4s$chromosome.x, .before = 1) %>% 
  dplyr::rename(chromosome = `plus_intersect_G4s$chromosome.x`) %>% 
  dplyr::rename(seq = `plus_intersect_G4s$seq`)

#save the gc fv

save(plus_intersect_gc, file = "plus_intersect_gc")

save.image("GC.RData")


# plus strand, absence of G4s

G4_exp_join_absence_plus <- anti_join(intersect_tibble, plus_intersect_G4s, by = c("seq" = "seq"))

g2_1 <- matrix(data = NA, nrow = dim(G4_exp_join_absence_plus)[1], ncol = 110)


for (j in 1:dim(G4_exp_join_absence_plus)[1]) {
  for (i in 1:110) {
    g2_1[j,i] <- GC(s2c(G4_exp_join_absence_plus$seq[j])[(chunks[i]):(chunks[i]+window_size-1)])
  }
  g2_1
}

# replace NAs

g2_1[is.na(g2_1)] <- 0


# assign to the right chromosome/sequences

plus_intersect_absence_gc <- tibble(as.data.table(g2_1)) %>% 
  add_column(G4_exp_join_absence_plus$seq, .before = 1) %>% 
  add_column(G4_exp_join_absence_plus$names, .before = 1) %>% 
  dplyr::rename(chromosome = `G4_exp_join_absence_plus$names`) %>% 
  dplyr::rename(seq = `G4_exp_join_absence_plus$seq`)

save(plus_intersect_absence_gc, file = "plus_intersect_absence_gc.RData")
