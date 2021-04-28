library(keras)
library(viridis)
library(magick)
library(tensorflow)
library(abind)

k_clear_session()
remove(model)

K <- backend()

use_compat(version = "v1")
setwd("/home/cjls4/ML/")

model <- load_model_hdf5("ML_2_model", 
                         compile = F)

graph <- tf$get_default_graph()

image <- get_file("elephant.jpg", "https://goo.gl/zCTWXW") %>% 
  image_load(target_size = c(224, 224)) %>% 
  image_to_array() %>% 
  array_reshape(dim = c(1, 224, 224, 3)) %>% 
  imagenet_preprocess_input()


dim(x_train) # 43826 2750 19

promoter <- x_train[200,,]

#x_train <- array_reshape(x_train, c(nrow(x_train), img_rows, img_cols, 1))

promoter <- array_reshape(promoter, c(1, 2750, 19))

tf$Graph$as_default(graph)

preds <- model %>% predict(promoter)
#dim num [1, 1:2], represents the 2 classes

#imagenet_decode_predictions(preds, top = 2)[[1]]

which.max(preds[1,])
#[1] 1, aka 0.645>0.355

tf$disable_eager_execution()

promoter_output <- model$output[, 1]
#Tensor("strided_slice:0", shape=(None,), dtype=float32)

last_conv_layer <- model %>% get_layer("conv1d")
#<tensorflow.python.keras.layers.convolutional.Conv1D>

grads <- K$gradients(promoter_output, last_conv_layer$output)[[1]]
# Tensor("gradients/average_pooling1d/ExpandDims_grad/Reshape:0", shape=(None, 2750, 1500), dtype=float32)

#something is wrong with the next line
pooled_grads <- K$mean(grads, axis = c(0L, 1L))
# Tensor("Mean:0", shape=(1500,), dtype=float32)

iterate <- K$`function`(list(model$input), 
                        list(pooled_grads, last_conv_layer$output[1,,])) 

# some of my pooled values are negative???
# well, a lot of them are
# uh.

c(pooled_grads_value, conv_layer_output_value) %<-% iterate(list(promoter))
# dim pooled_grads_value 1500
# dim conv_layer_output_value 2750 1500


model$input
last_conv_layer$output

for (i in 1:1500) {
  conv_layer_output_value[,i] <- 
    conv_layer_output_value[,i] * pooled_grads_value[[i]] 
}

#heatmap <- apply(conv_layer_output_value, c(1,2), mean)
heatmap <- apply(conv_layer_output_value, 1, mean)

hist(heatmap)

heatmap <- heatmap*10

heatmap <- pmax(heatmap, 0) 
heatmap <- heatmap / max(heatmap)

write_heatmap <- function(heatmap, filename, width = 1, height = 2750,
                          bg = "white", col = terrain.colors(12)) {
  png(filename, width = width, height = height, bg = bg, type = "cairo")
  op = par(mar = c(0,0,0,0))
  on.exit({par(op); dev.off()}, add = TRUE)
  rotate <- function(x) t(apply(x, 2, rev))
  image(rotate(heatmap), axes = FALSE, asp = 1, col = col)
}


write_heatmap_2 <- function(heatmap, filename, width = 250, height = 2750,
                          bg = "white", col = terrain.colors(12)) {
  png(filename, width = width, height = height, bg = bg, type = "cairo")
  op = par(mar = c(0,0,0,0))
  on.exit({par(op); dev.off()}, add = TRUE)
  #rotate <- function(x) t(apply(x, 2, rev))
  image(t(as.matrix(heatmap)), axes = FALSE, col = col)
}

write_heatmap_2(heatmap, "second_try_heatmap.png") 

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

write_heatmap_3 <- function(heatmap, filename, width = 1000, height = 2750,
                            bg = "white", col = terrain.colors(12)) {
  png(filename, width = width, height = height, bg = bg, type = "cairo")
  op = par(mar = c(0,0,0,0))
  on.exit({par(op); dev.off()}, add = TRUE)
  #rotate <- function(x) t(apply(x, 2, rev))
  image(t(as.matrix(heatmap)), axes = FALSE, col = col)
}

write_heatmap_3(heatmap3, "third_try_heatmap.png") 

write_heatmap_3(heatmap2, "4th_try_hm.png")

load("x_train_min_p")

bruh <- x_train_min_p[,,200]

bruh_hm <- cbind(heatmap, bruh)

bruh_hm[is.na(bruh_hm)] <- 0

#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#apply Min-Max normalization to first four columns in iris dataset
bruh_norm <- as.data.frame(lapply( as.data.frame(bruh_hm), min_max_norm))

bruh_norm[is.na(bruh_norm)] <- 0
bruh_norm <- cbind(bruh_norm, t(specific))

setwd("/home/cjls4/ML/")

write_heatmap_3(bruh_norm, "5th_try.png")

bruh_check <- cbind(rev(t(heatmap)), bruh)
bruh_check[is.na(bruh_check)] <- 0
bruh_check <- as.data.frame(lapply(as.data.frame(bruh_check), min_max_norm))

bruh_check <- cbind(bruh_check, t(specific))

write_heatmap_3(bruh_check, "6th_try.png")

# can try and do correlations across different promoters and the 

setwd("/home/cjls4/feature_vectors/")
