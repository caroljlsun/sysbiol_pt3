library(corrplot)
library(tidyverse)
library(stats)
library(useful)

setwd("/home/cjls4/feature_vectors/")

load("x_neg.RData")

##### generic explorations of the data #####
dim(x_total)

#average a content
average_a <- x_neg[,,1]

average_a <- apply(average_a, 2, mean)

avg_a <- cbind(tibble(average_a), tibble(c(1:2750)))

ggplot(avg_a, mapping = aes(x = `c(1:2750)`, 
                            y = average_a,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average A content")+
  ggtitle("Average A content, G4 negative")

#average c content
average_c <- x_neg[,,2]

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
  ggtitle("Average C content, G4 negative")


#average g content
average_g <- x_neg[,,3]
average_g <- apply(average_g, 2, mean)

avg_g <- cbind(tibble(average_g), tibble(c(1:2750)))

ggplot(avg_g, mapping = aes(x = `c(1:2750)`, 
                            y = average_g,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average G content")+
  ggtitle("Average G content, G4 negative")



#average t content
average_t <- x_neg[,,4]
average_t <- apply(average_t, 2, mean)

avg_t <- cbind(tibble(average_t), tibble(c(1:2750)))

ggplot(avg_t, mapping = aes(x = `c(1:2750)`, 
                            y = average_t,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average T content")+
  ggtitle("Average T content, G4 negative")

# average histone marks

#h3k27ac_b

average_h3k27ac_b <- x_neg[,,9]
average_h3k27ac_b <- apply(average_h3k27ac_b, 2, mean)

avg_h3k27ac_b <- cbind(tibble(average_h3k27ac_b), tibble(c(1:2750)))

ggplot(avg_h3k27ac_b, mapping = aes(x = `c(1:2750)`, 
                                    y = average_h3k27ac_b,
                                    group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.25, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average H3K27ac_b score")+
  ggtitle("Average H3K27ac_b distribution, G4 negative")

# 

average_h3k27ac_t <- x_neg[,,10]
average_h3k27ac_t <- apply(average_h3k27ac_t, 2, mean)

avg_h3k27ac_t <- cbind(tibble(average_h3k27ac_t), tibble(c(1:2750)))

ggplot(avg_h3k27ac_t, mapping = aes(x = `c(1:2750)`, 
                                    y = average_h3k27ac_t,
                                    group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.25, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average H3K27ac_t score")+
  ggtitle("Average H3K27ac_t distribution, G4 negative")

#
average_h3k27me3_b <- x_neg[,,11]
average_h3k27me3_b <- apply(average_h3k27me3_b, 2, mean)

avg_h3k27me3_b <- cbind(tibble(average_h3k27me3_b), tibble(c(1:2750)))

ggplot(avg_h3k27me3_b, mapping = aes(x = `c(1:2750)`, 
                                     y = average_h3k27me3_b,
                                     group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.25, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average H3K27me3_b score")+
  ggtitle("Average H3K27me3_b distribution, G4 negative")

#h3k27me3_t

average_h3k27me3_t <- x_neg[,,12]
average_h3k27me3_t <- apply(average_h3k27me3_t, 2, mean)

avg_h3k27me3_t <- cbind(tibble(average_h3k27me3_t), tibble(c(1:2750)))

ggplot(avg_h3k27me3_t, mapping = aes(x = `c(1:2750)`, 
                                     y = average_h3k27me3_t,
                                     group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.25, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average H3K27me3_t score")+
  ggtitle("Average H3K27me3_t distribution, G4 negative")

# average r loop positions

average_r_loops <- x_neg[,,13]
average_r_loops <- apply(average_r_loops, 2, mean)

avg_r <- cbind(tibble(average_r_loops), tibble(c(1:2750)))

ggplot(avg_r, mapping = aes(x = `c(1:2750)`, 
                            y = average_t,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average R-loop score")+
  ggtitle("Average R-loop position, G4 negative")

# average methylation

average_methylation <- x_neg[,,18]
average_methylation <- apply(average_methylation, 2, mean)

avg_m <- cbind(tibble(average_methylation), tibble(c(1:2750)))

ggplot(avg_m, mapping = aes(x = `c(1:2750)`, 
                            y = average_methylation,
                            group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.3, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average Methylation score")+
  ggtitle("Average methylation distribution, G4 negative")

# average gc content across positions

average_gc <- x_neg[,,19]
average_gc <- apply(average_gc, 2, mean)

avg_gc <- cbind(tibble(average_gc), tibble(c(1:2750)))

ggplot(avg_gc, mapping = aes(x = `c(1:2750)`, 
                             y = average_gc,
                             group = 1))+
  geom_point(col = "grey", alpha = 0.5)+
  geom_smooth(method = "loess", span = 0.2, col = "red", fill = "cadetblue")+
  theme_bw()+
  xlab("Position")+
  ylab("Average GC score")+
  ggtitle("Average GC distribution, G4 negative")

