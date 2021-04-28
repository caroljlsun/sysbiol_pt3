
library(keras)
library(useful)
library(abind)
library(tensorflow)
library(plyr)
library(data.table)

# global options

setwd("/home/cjls4/feature_vectors/")

# loading data

load("x_train_min_p")
load("x_train_min_n")
load("x_train_plus_p")
load("x_train_plus_n")

load("y_train_min_p")
load("y_train_min_n")
load("y_train_plus_p")
load("y_train_plus_n")

load("x_test_min_p")
load("x_test_min_n")
load("x_test_plus_p")
load("x_test_plus_n")

load("y_test_min_p")
load("y_test_min_n")
load("y_test_plus_p")
load("y_test_plus_n")


##### temporary fix to the NA's problem #####

g4_positive <- abind(x_train_min_p, x_train_plus_p)

g4_positive[is.na(g4_positive)] <- 0

g4_negative <- abind(x_train_min_n, x_train_plus_n)

g4_negative[is.na(g4_negative)] <- 0

dim(g4_positive) #2750 19 10401 
##### turn 3D array into 2D matrix #####

collapsed_g4_positive <- apply(g4_positive, c(2,3), mean)

collapsed_g4_positive <- as.data.table(t(collapsed_g4_positive))

colnames(collapsed_g4_positive) <- c("A", "C", "G", "T", "bcellf", "bcellm", "tcellf", "tcellm", "H3K27ac_B",
                                     "H3K27ac_T", "H3K27me3_B", "H3K27me3_T", "R_loops", "MGW", "ProT", 
                                     "Roll", "HelT", "Methylation", "GC")

collapsed_g4_positive$G4 <- 1


#g4 negatives

collapsed_g4_negative <- apply(g4_negative, c(2,3), mean)

collapsed_g4_negative <-  as.data.table(t(collapsed_g4_negative))

colnames(collapsed_g4_negative) <- c("A", "C", "G", "T", "bcellf", "bcellm", "tcellf", "tcellm", "H3K27ac_B",
                                     "H3K27ac_T", "H3K27me3_B", "H3K27me3_T", "R_loops", "MGW", "ProT", 
                                     "Roll", "HelT", "Methylation", "GC")

collapsed_g4_negative$G4 <- 0


#####  now for the test versions #####


test_g4_positive <- abind(x_test_min_p, x_test_plus_p)

test_g4_positive[is.na(test_g4_positive)] <- 0

test_g4_negative <- abind(x_test_min_n, x_test_plus_n)

test_g4_negative[is.na(test_g4_negative)] <- 0

dim(test_g4_positive) #2750 19 10401 
##### turn 3D array into 2D matrix #####

test_collapsed_g4_positive <- apply(test_g4_positive, c(2,3), mean)

test_collapsed_g4_positive <- as.data.table(t(test_collapsed_g4_positive))

colnames(test_collapsed_g4_positive) <- c("A", "C", "G", "T", "bcellf", "bcellm", "tcellf", "tcellm",
                                          "H3K27ac_B", "H3K27ac_T", "H3K27me3_B", "H3K27me3_T", "R_loops",
                                          "MGW", "ProT", "Roll", "HelT", "Methylation", "GC")

test_collapsed_g4_positive$G4 <- 1


#g4 negatives

test_collapsed_g4_negative <- apply(test_g4_negative, c(2,3), mean)

test_collapsed_g4_negative <-  as.data.table(t(test_collapsed_g4_negative))

colnames(test_collapsed_g4_negative) <- c("A", "C", "G", "T", "bcellf", "bcellm", "tcellf", "tcellm",
                                          "H3K27ac_B", "H3K27ac_T", "H3K27me3_B", "H3K27me3_T", "R_loops",
                                          "MGW", "ProT", "Roll", "HelT", "Methylation", "GC")

test_collapsed_g4_negative$G4 <- 0


#save objects

save(collapsed_g4_positive, file = "collapsed_g4_positive.RData")

save(collapsed_g4_negative, file = "collapsed_g4_negative.RData")

save(test_collapsed_g4_positive, file = "test_collapsed_g4_positive.RData")

save(test_collapsed_g4_negative, file = "test_collapsed_g4_negative.RData")
