library(keras)
K <- backend()

model <- application_vgg16(weights = "imagenet") 
graph <- tf$get_default_graph()
image <- get_file("elephant.jpg", "https://goo.gl/zCTWXW") %>% 
  image_load(target_size = c(224, 224)) %>% 
  image_to_array() %>% 
  array_reshape(dim = c(1, 224, 224, 3)) %>% 
  imagenet_preprocess_input()

tf$Graph$as_default(graph)
preds <- model %>% predict(image)

#dim preds [1, 1:1000]


imagenet_decode_predictions(preds, top = 3)[[1]]

#  class_name class_description       score
# 1  n02504458  African_elephant 0.909421027
# 2  n01871265            tusker 0.086182885
# 3  n02504013   Indian_elephant 0.004354583


uhhh <-  get_file("elephant.jpg", "https://goo.gl/zCTWXW") %>% 
  image_load(target_size = c(224, 224))

im <- image[1,,,1]

image(x=1:224, y=1:224, im, col=gray((0:255)/255))

image_write(uhhh, path = "check.png")

# class_name class_description      score
# 1  n02504458  African_elephant 0.78988522
# 2  n01871265            tusker 0.19872670
# 3  n02504013   Indian_elephant 0.01114247

which.max(preds[1,])

# [1] 387
# 
# tf.compat.v1.disable_eager_execution()
# 
# ?tf_extract_opts
# 
# ?tf$GradientTape$gradient
# tf$GradientTape$watch (african_elephant_output, last_conv_layer$output)
# Get gradient of the winner class w.r.t. the output of the (last) conv. layer
# with tf.GradientTape() as gtape:
#   conv_output, predictions = heatmap_model(img_tensor)
# loss = predictions[:, np.argmax(predictions[0])]
# grads = gtape.gradient(loss, conv_output)
# pooled_grads = K.mean(grads, axis=(0, 1, 2))
# 
# heatmap = tf.reduce_mean(tf.multiply(pooled_grads, conv_output), axis=-1)
# heatmap = np.maximum(heatmap, 0)
# max_heat = np.max(heatmap)
# if max_heat == 0:
#   max_heat = 1e-10
# heatmap /= max_heat
# 
# print(heatmap.shape)

tf$disable_eager_execution()

african_elephant_output <- model$output[, 387]
# Tensor("strided_slice:0", shape=(None,), dtype=float32)

last_conv_layer <- model %>% get_layer("block5_conv3")
#<tensorflow.python.keras.layers.convolutional.Conv2D>


grads <- K$gradients(african_elephant_output, last_conv_layer$output)[[1]]
#Tensor("gradients/block5_pool/MaxPool_grad/MaxPoolGrad:0", shape=(None, 14, 14, 512), dtype=float32)

pooled_grads <- K$mean(grads, axis = c(0L, 1L, 2L))
#Tensor("Mean:0", shape=(512,), dtype=float32)

iterate <- K$`function`(list(model$input), 
                        list(pooled_grads, last_conv_layer$output[1,,,])) 


c(pooled_grads_value, conv_layer_output_value) %<-% iterate(list(image))
# dim pooled_grads_value num [1:512(1d)], or [1] 512
# dim conv_layer_output_value num [1:14, 1:14, 1:512]

model$input
# Tensor("input_1:0", shape=(None, 224, 224, 3), dtype=float32)


for (i in 1:512) {
  conv_layer_output_value[,,i] <- 
    conv_layer_output_value[,,i] * pooled_grads_value[[i]] 
}
# 14 14 512


heatmap <- apply(conv_layer_output_value, c(1,2), mean)
# X : an array
# MARGIN: a vector giving the subscripts which the function will be applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns.
# FUN :  the function
# dim num [1:14, 1:14]
# the literal final size of the heatmap

heatmap <- pmax(heatmap, 0) 
# use pmin(5:1, pi) to get what it does

heatmap <- heatmap / max(heatmap)

write_heatmap <- function(heatmap, filename, width = 224, height = 224,
                          bg = "white", col = terrain.colors(12)) {
  png(filename, width = width, height = height, bg = bg, type = "cairo")
  op = par(mar = c(0,0,0,0))
  on.exit({par(op); dev.off()}, add = TRUE)
  rotate <- function(x) t(apply(x, 2, rev))
  image(rotate(heatmap), axes = FALSE, asp = 1, col = col)
}

write_heatmap(heatmap, "elephant_heatmap.png") 

library(magick) 
library(viridis) 


image <- image_read("elephant_heatmap.png") 
info <- image_info(image) 
geometry <- sprintf("%dx%d!", info$width, info$height) 

pal <- col2rgb(viridis(20), alpha = TRUE)
alpha <- floor(seq(0, 255, length = ncol(pal))) 
pal_col <- rgb(t(pal), alpha = alpha, maxColorValue = 255)
write_heatmap(heatmap, "elephant_overlay.png", 
              width = 14, height = 14, bg = NA, col = pal_col) 

image_read("elephant_overlay.png") %>% 
  image_resize(geometry, filter = "quadratic") %>% 
  image_composite(image, operator = "blend", compose_args = "20") %>%
  plot() 

