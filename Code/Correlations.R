library(corrplot)
library(tidyverse)
library(stats)
library(useful)

setwd("/home/cjls4/feature_vectors/")

load("collapsed_g4_positive.RData")
load("collapsed_g4_negative.RData")
load("test_collapsed_g4_positive.RData")
load("test_collapsed_g4_negative.RData")

load("minus_intersect_G4s.RData")
load("plus_intersect_G4s.RData")
#merge the +/- training sets

train_g4 <- rbind(collapsed_g4_positive, collapsed_g4_negative)


#merge the +/- test sets

test_g4 <-rbind(test_collapsed_g4_positive, test_collapsed_g4_negative)


total_g4 <- rbind(train_g4, test_g4)

#GC_and_R_loops <- cor(x_train[,,19,1], x_train[,,13,1])

collapsed_g4_positive$G4

lol <- cor(collapsed_g4_positive[,-20], collapsed_g4_positive[,-20])
corrplot::corrplot(lol, tl.col = "black")

lol2 <- cor(collapsed_g4_negative[,-20], collapsed_g4_negative[,-20])
corrplot::corrplot(lol2, tl.col = "black")

#barely any differences lol


# loading data

load("min_g4_positive_names_train.RData")
load("min_g4_negative_names_train.RData")
load("plus_g4_positive_names_train.RData")
load("plus_g4_negative_names_train.RData")

load("min_g4_positive_names_test.RData")
load("min_g4_negative_names_test.RData")
load("plus_g4_positive_names_test.RData")
load("plus_g4_negative_names_test.RData")

load("x_train_min_p")
load("x_train_min_n")
load("x_train_plus_p")
load("x_train_plus_n")

load("x_test_min_p")
load("x_test_min_n")
load("x_test_plus_p")
load("x_test_plus_n")

load("collapsed_g4_positive.RData")
load("collapsed_g4_negative.RData")
load("test_collapsed_g4_positive.RData")
load("test_collapsed_g4_negative.RData")


#merge the +/- training sets

train_g4 <- rbind(collapsed_g4_positive, collapsed_g4_negative)

remove(collapsed_g4_positive)
remove(collapsed_g4_negative)

#merge the +/- test sets

test_g4 <-rbind(test_collapsed_g4_positive, test_collapsed_g4_negative)

remove(test_collapsed_g4_positive)
remove(test_collapsed_g4_negative)

total_g4 <- rbind(train_g4, test_g4)


x_train <- abind(x_train_min_p, x_train_min_n, x_train_plus_p, x_train_plus_n) 
x_test <- abind(x_test_min_p, x_test_min_n, x_test_plus_p, x_test_plus_n)

dim(x_train) #2750 19 43826

remove(x_train_min_p)
remove(x_train_min_n)
remove(x_train_plus_p)
remove(x_train_plus_n)


x_names <- rbind(min_g4_positive_names_train, min_g4_negative_names_train,
                 plus_g4_positive_names_train, plus_g4_negative_names_train,
                 min_g4_positive_names_test, min_g4_negative_names_test,
                 plus_g4_positive_names_test, plus_g4_negative_names_test)

##### temporary fix to the NA's problem #####

x_train[is.na(x_train)] <- 0
x_test[is.na(x_test)] <- 0

x_train <- aperm(x_train, c(3,1,2))
x_test <- aperm(x_test, c(3,1,2))
#normalise

x_train <- keras::normalize(x_train, axis = -1, order = 2)
x_test <- keras::normalize(x_test, axis = -1, order = 2)

#combine x_train with x_test

x_total <- abind(x_train, x_test, along = 1)

x_pos <- x_total[1:500,,]
x_neg <- x_total[6000:6500,,]

setwd("/home/cjls4/feature_vectors/")
save(x_pos, file = "x_pos.RData")
save(x_neg, file = "x_neg.RData")

load("x_pos.RData")

##### generic explorations of the data #####
dim(x_total)

#average a content
average_a <- x_pos[,,1]

average_a <- apply(average_a, 2, mean)

plot(average_a, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_a <- cbind(tibble(average_a), tibble(c(1:2750)))

ggplot(avg_a, mapping = aes(x = `c(1:2750)`, 
                            y = average_a,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average A content")+
  ggtitle("Average A content, G4 positive")

#average c content
average_c <- x_pos[,,2]

average_c <- apply(average_c, 2, mean)

avg_c <- cbind(tibble(average_c), tibble(c(1:2750)))

ggplot(avg_c, mapping = aes(x = `c(1:2750)`, 
                            y = average_c,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average C content")+
  ggtitle("Average C content, G4 positive")


#average g content
average_g <- x_pos[,,3]
average_g <- apply(average_g, 2, mean)

plot(average_g, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_g <- cbind(tibble(average_g), tibble(c(1:2750)))

ggplot(avg_g, mapping = aes(x = `c(1:2750)`, 
                            y = average_g,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average G content")+
  ggtitle("Average G content, G4 positive")



#average t content
average_t <- x_pos[,,4]
average_t <- apply(average_t, 2, mean)

plot(average_t, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_t <- cbind(tibble(average_t), tibble(c(1:2750)))

ggplot(avg_t, mapping = aes(x = `c(1:2750)`, 
                            y = average_t,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average T content")+
  ggtitle("Average T content, G4 positive")


# average g4 positions 

average_g4s <- apply(plus_intersect_G4s[,3:2752], 2, mean)

plot(average_g4s, type = "l", xlab = "Position", ylab = "G4 score", main = "Average G4 distribution")


avg_g4s <- cbind(tibble(average_g4s), tibble(c(2750:1)))

ggplot(avg_g4s, mapping = aes(x = `c(2750:1)`, 
                            y = average_g4s,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.09, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average G4 score")+
  ggtitle("Average G4 position, G4 positive")

# average histone marks

#h3k27ac_b

average_h3k27ac_b <- x_pos[,,9]
average_h3k27ac_b <- apply(average_h3k27ac_b, 2, mean)

plot(average_h3k27ac_b, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_h3k27ac_b <- cbind(tibble(average_h3k27ac_b), tibble(c(1:2750)))

ggplot(avg_h3k27ac_b, mapping = aes(x = `c(1:2750)`, 
                            y = average_h3k27ac_b,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.25, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average H3K27ac_b score")+
  ggtitle("Average H3K27ac_b distribution, G4 positive")

# 

average_h3k27ac_t <- x_pos[,,10]
average_h3k27ac_t <- apply(average_h3k27ac_t, 2, mean)

plot(average_h3k27ac_t, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_h3k27ac_t <- cbind(tibble(average_h3k27ac_t), tibble(c(1:2750)))

ggplot(avg_h3k27ac_t, mapping = aes(x = `c(1:2750)`, 
                                    y = average_h3k27ac_t,
                                    group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.25, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average H3K27ac_t score")+
  ggtitle("Average H3K27ac_t distribution, G4 positive")

#
average_h3k27me3_b <- x_pos[,,11]
average_h3k27me3_b <- apply(average_h3k27me3_b, 2, mean)

plot(average_h3k27me3_b, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_h3k27me3_b <- cbind(tibble(average_h3k27me3_b), tibble(c(1:2750)))

ggplot(avg_h3k27me3_b, mapping = aes(x = `c(1:2750)`, 
                                    y = average_h3k27me3_b,
                                    group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.25, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average H3K27me3_b score")+
  ggtitle("Average H3K27me3_b distribution, G4 positive")

#h3k27me3_t

average_h3k27me3_t <- x_pos[,,12]
average_h3k27me3_t <- apply(average_h3k27me3_t, 2, mean)

plot(average_h3k27me3_t, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_h3k27me3_t <- cbind(tibble(average_h3k27me3_t), tibble(c(1:2750)))

ggplot(avg_h3k27me3_t, mapping = aes(x = `c(1:2750)`, 
                                     y = average_h3k27me3_t,
                                     group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.25, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average H3K27me3_t score")+
  ggtitle("Average H3K27me3_t distribution, G4 positive")

# average r loop positions

average_r_loops <- x_pos[,,13]
average_r_loops <- apply(average_r_loops, 2, mean)

plot(average_r_loops, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_r <- cbind(tibble(average_r_loops), tibble(c(1:2750)))

ggplot(avg_r, mapping = aes(x = `c(1:2750)`, 
                            y = average_t,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average R-loop score")+
  ggtitle("Average R-loop position, G4 positive")

# average methylation

average_methylation <- x_pos[,,18]
average_methylation <- apply(average_methylation, 2, mean)

plot(average_methylation, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_m <- cbind(tibble(average_methylation), tibble(c(1:2750)))

ggplot(avg_m, mapping = aes(x = `c(1:2750)`, 
                            y = average_methylation,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average Methylation score")+
  ggtitle("Average methylation distribution, G4 positive")

# average gc content across positions

average_gc <- x_pos[,,19]
average_gc <- apply(average_gc, 2, mean)

plot(average_gc, type = "p", xlab = "Position", ylab = "R-loops score", main = "Average R-loop distribution")

avg_gc <- cbind(tibble(average_gc), tibble(c(1:2750)))

ggplot(avg_gc, mapping = aes(x = `c(1:2750)`, 
                            y = average_gc,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.2, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average GC score")+
  ggtitle("Average GC distribution, G4 positive")

