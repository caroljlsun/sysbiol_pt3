library(keras)
library(viridis)
library(magick)
library(tensorflow)
library(abind)

k_clear_session()
remove(model)

use_compat(version = "v1")
setwd("/home/cjls4/ML/")

# model <- load_model_hdf5("ML_4_model", 
#                          compile = F)

model <- load_model_tf("ML_4_model", compile = F)

K <- backend()
graph <- tf$get_default_graph()


dim(x_total) # 46528 2750 19

#promoter <- x_total

x_total <- array_reshape(x_total, c(nrow(x_total), 2750, 19))

#x_total <- array_reshape(x_total, c(1, 2750, 19))

tf$Graph$as_default(graph)

all_preds <- model %>% predict(x_train)
#dim num [1, 1:2], represents the 2 classes

#imagenet_decode_predictions(preds, top = 2)[[1]]

which.max(all_preds[1,])

#[1] 12672, aka 0.9998881

tf$disable_eager_execution()

promoter_output <- model$output[,1]
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

pgv1 <- array(data = NA, c(1200,1500))
clov1 <- array(data = NA, c(1200,2750,1500))

#the following numbers were chosen to be 500 promoters from positive G4s, minus strand, plus strand, train and
#test set

for (i in c(1:500)) {
  
  test <- x_total[i,,]
  test <- array_reshape(test, c(1, 2750, 19))
  c(pgv1[i,], clov1[i,,]) %<-% iterate(list(test))
  
}

# , 21914:22413, 43827:43926, 45178:45277

#c(pooled_grads_value, conv_layer_output_value) %<-% iterate(list(test))


# dim pooled_grads_value 1500
# dim conv_layer_output_value 2750 1500


model$input
last_conv_layer$output

for (j in 1:500) {
  for (i in 1:1500) {
    clov1[j,,i] <- clov1[j,,i]*pgv1[j,i]
  }
}


for (i in 1:1500) {
  conv_layer_output_value[,i] <- 
    conv_layer_output_value[,i] * pooled_grads_value[[i]] 
}

heatmaps_500 <- apply(clov1[1:500,,], 2, mean)

#save(heatmaps_500, file = "heatmaps_500.RData")
load("/home/cjls4/feature_vectors/heatmaps_500.RData")
plot(heatmaps_500)

heatmaps_500 <- pmax(heatmaps_500, 0) 
heatmaps_500 <- heatmaps_500 / max(heatmaps_500)
hist(heatmaps_500)


write_heatmap_2 <- function(heatmap, filename, width = 250, height = 2750,
                            bg = "white", col = viridis(12, begin = 0, end = 0.9)) {
  png(filename, width = width, height = height, bg = bg, type = "cairo")
  op = par(mar = c(0,0,0,0))
  on.exit({par(op); dev.off()}, add = TRUE)
  #rotate <- function(x) t(apply(x, 2, rev))
  image(t(as.matrix(heatmap)), axes = FALSE, col = col)
}


write_heatmap_2(heatmaps_500, "tryouts.png") 

#### ignore for now ####

library(tidyverse)

x_train_min_p

load("minus_intersect_G4s.RData")
load("expression_tibble.RData")

G4_exp_join_min <- inner_join(intersect_tibble, minus_intersect_G4s, by = c("seq" = "seq"))
full_minimal <- G4_exp_join_min

p_chr1s <- str_detect(full_minimal$names, "chr1_")

p_chr1s_opp <- str_detect(full_minimal$names, "chr1_", negate = TRUE)

p_chr1_g4s <- G4_exp_join_min[,9:2758]

p_chr1_g4s_avg <- apply(p_chr1_g4s, 2, mean)

plot(p_chr1_g4s_avg)

p_chr1_g4s_avg <- pmax(p_chr1_g4s_avg, 0) 
p_chr1_g4s_avg <- p_chr1_g4s_avg / max(p_chr1_g4s_avg)
hist(p_chr1_g4s_avg)


write_heatmap_2(p_chr1_g4s_avg, "average_g4s.png") 





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
