# same as before but for the one hot encoding only

library(keras)
library(viridis)
library(magick)
library(tensorflow)
library(abind)

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

x_train <- abind(x_train_min_p, x_train_min_n, x_train_plus_p, x_train_plus_n)
y_train <- c(y_train_min_p, y_train_min_n, y_train_plus_p, y_train_plus_n)

remove(x_train_min_p)
remove(x_train_min_n)
remove(x_train_plus_p)
remove(x_train_plus_n)

x_test <- abind(x_test_min_p, x_test_min_n, x_test_plus_p, x_test_plus_n)
y_test <- c(y_test_min_p, y_test_min_n, y_test_plus_p, y_test_plus_n)

##### one_hot excluded #####

x_train_min <- x_train[,1:4,]
y_train_min <- y_train

x_test_min <- x_test[,1:4,]
y_test_min <- y_test

remove(x_train)
remove(y_train)

##### temporary fix to the NA's problem #####

x_train_min[is.na(x_train_min)] <- 0
x_test_min[is.na(x_test_min)] <- 0

y_train_min <- to_categorical(y_train_min, 2)
y_test_min <- to_categorical(y_test_min, 2)

#dimensions

dim(x_train_min) # 2750 4 43826 (rows, col, slices)
dim(y_train_min) # 43826 

#melt some arrays

x_train_min <- aperm(x_train_min, c(3,1,2))
x_test_min <- aperm(x_test_min, c(3,1,2))


# check dimensions again

dim(x_train_min) # 43826 2750 4 (slices, row, col)
dim(y_train_min) # 43826 2  (slices)

##### normalise #####

x_train_min <- keras::normalize(x_train_min, axis = -1, order = 2)
x_test_min <- keras::normalize(x_test_min, axis = -1, order = 2)

##### extract weights etc #####

k_clear_session()
remove(model)

use_compat(version = "v1")
setwd("/home/cjls4/ML/")

# model <- load_model_hdf5("ML_5_87AUC_model/", 
#                          compile = F)

model <- load_model_tf("ML_5_model", compile = F)

K <- backend()
graph <- tf$get_default_graph()



dim(x_train_min) # 46528 2750 15

#promoter <- x_train_min

#x_train_min <- array_reshape(x_train_min, c(nrow(x_train_min), 2750, 15))

#x_train_min <- array_reshape(x_train_min, c(1, 2750, 19))

tf$Graph$as_default(graph)

all_preds <- model %>% predict(x_train_min)
#dim num [1, 1:2], represents the 2 classes

#This is needed to get gradients, strangely difficult to find function
tf$disable_eager_execution()

promoter_output <- model$output[,1]
#Tensor("strided_slice:0", shape=(None,), dtype=float32)

last_conv_layer <- model %>% get_layer(index = 1)
#<tensorflow.python.keras.layers.convolutional.Conv1D>

grads <- K$gradients(promoter_output, last_conv_layer$output)[[1]]
# Tensor("gradients/average_pooling1d/ExpandDims_grad/Reshape:0", shape=(None, 2750, 1500), dtype=float32)

pooled_grads <- K$mean(grads, axis = c(0L, 1L))
# Tensor("Mean:0", shape=(1500,), dtype=float32)

iterate <- K$`function`(list(model$input), 
                        list(pooled_grads, last_conv_layer$output[1,,])) 

pgv1 <- array(data = NA, c(500,1500))
clov1 <- array(data = NA, c(500,2750,1500))

#the following numbers were chosen to be 500 promoters from positive G4s, minus strand. More than this takes too long

for (i in c(1:500)) {
  
  test <- x_train_min[i,,]
  test <- array_reshape(test, c(1, 2750, 4))
  c(pgv1[i,], clov1[i,,]) %<-% iterate(list(test))
  
}


# multiply the clov by pgv

for (j in 1:500) {
  for (i in 1:1500) {
    clov1[j,,i] <- clov1[j,,i]*pgv1[j,i]
  }
}



heatmaps_model_5 <- apply(clov1[1:500,,], 2, mean)

plot(heatmaps_model_5)

save(heatmaps_model_5, file = "heatmaps_model_5.RData")
load("ML/heatmaps_model_5.RData")


heatmaps_model_5 <- pmax(heatmaps_model_5, 0) 
heatmaps_model_5 <- heatmaps_model_5 / max(heatmaps_model_5)
hist(heatmaps_model_5)


write_heatmap_2 <- function(heatmap, filename, width = 250, height = 2750,
                            bg = "white", col = viridis(12, begin = 0, end = 0.9)) {
  png(filename, width = width, height = height, bg = bg, type = "cairo")
  op = par(mar = c(0,0,0,0))
  on.exit({par(op); dev.off()}, add = TRUE)
  #rotate <- function(x) t(apply(x, 2, rev))
  image(t(as.matrix(heatmap)), axes = FALSE, col = col)
}


write_heatmap_2(heatmaps_model_5, "tryouts_feature_masked_one_hot.png") 


##### ignore for now #####

library(tidyverse)

x_train_min_p

load("minus_intersect_G4s.RData")
load("expression_tibble.RData")

G4_exp_join_min <- inner_join(intersect_tibble, minus_intersect_G4s, by = c("seq" = "seq"))


specific <- p_chr1_G4s[200,]

max(heatmap)

specific[specific == 1] <- max(heatmap)

heatmap2 <- cbind(heatmap, t(specific))

heatmap3 <- cbind(heatmap, rev(t(specific)))



setwd("/home/cjls4/feature_vectors/")
