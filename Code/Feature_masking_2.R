
#remove which features?

# R loops
# GC content
# Purely one hot encoding?

# libraries

library(keras)
library(useful)
library(abind)
library(tensorflow)

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

# ml

batch_size <- 5
num_classes <- 2
epochs <- 15

img_rows <- 2750
img_cols <- 15

#

x_train <- abind(x_train_min_p, x_train_min_n, x_train_plus_p, x_train_plus_n)
y_train <- c(y_train_min_p, y_train_min_n, y_train_plus_p, y_train_plus_n)

remove(x_train_min_p)
remove(x_train_min_n)
remove(x_train_plus_p)
remove(x_train_plus_n)

x_test <- abind(x_test_min_p, x_test_min_n, x_test_plus_p, x_test_plus_n)
y_test <- c(y_test_min_p, y_test_min_n, y_test_plus_p, y_test_plus_n)

##### one_hot excluded #####

x_train_min <- x_train[,5:19,]
y_train_min <- y_train

x_test_min <- x_test[,5:19,]
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


input_shape <- c(img_rows, img_cols)


# the deepG4 paper uses a ff nn
# 1st layer 1D convolutional layer, 900 kernels, kernel size 20bp
# 2nd layer local average pooling, pool size 10bp
# 3rd layer global max pooling
# 4th layer drop out for regularisation
# 5th layer dense layer 
# 6th layer activation layer, activation sigmoid

#Hyperparameters including the number of kernels (900), kernel size (20bp), kernel activation (relu),
#pool size (10bp), drop-out (0%), epoch number (20), number of neurons in the dense layer (100) and the optimizer choice (rmsprop) were found by fine-tuning with a random grid.

model <- keras_model_sequential() %>% 
  layer_conv_1d(filters = 1500, kernel_size = 20, activation = 'relu',
                input_shape = input_shape, padding = "same") %>% 
  layer_average_pooling_1d(pool_size = 10) %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0) %>% 
  layer_dense(units = 100) %>% 
  layer_dense(units = num_classes, activation = 'sigmoid')

# Compile model
model %>% compile(
  loss = loss_binary_crossentropy,
  optimizer = optimizer_adam(),
  metrics = tf$keras$metrics$AUC()
)

# Train model
# model %>% fit(
#   x_train_min, y_train_min,
#   batch_size = batch_size,
#   epochs = epochs,
#   validation_split = 0.2
# )

cnn <- fit(model,
            x_train_min, y_train_min,
            batch_size = batch_size,
            epochs = epochs,
            validation_split = 0.2)

model <- load_model_hdf5("ML_6_81AUC_model/", compile = F)

scores <- model %>% evaluate(
  x_test_min, y_test_min, verbose = 1
)

# Output metrics
cat('Test loss:', scores[[1]], '\n')
cat('Test accuracy:', scores[[2]], '\n')

#Test loss: 0.5299206 
#Test accuracy: 0.8101588  

hm <- model$history

plot(
  cnn,
  metrics = NULL,
  method = c("auto", "ggplot2", "base"),
  smooth = getOption("keras.plot.history.smooth", TRUE),
  theme_bw = getOption("keras.plot.history.theme_bw", FALSE)
)



#save model
setwd("/home/cjls4/ML/")
save_model_hdf5(object = model, file = "ML_6_model")
save_model_tf(model, filepath = "ML_6_81AUC_model", include_optimizer = T)

#use this to delete old models

#remove(model)

#remove(scores)

#k_clear_session()
