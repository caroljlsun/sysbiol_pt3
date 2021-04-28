library(keras)
library(viridis)
library(tensorflow)
library(caret)
library(DESeq2)

k_clear_session()
remove(model)

#global options

use_compat(version = "v1")
setwd("/home/cjls4/ML/")

model <- load_model_hdf5("ML_4_model", 
                         compile = F)

graph <- tf$get_default_graph()

tf$Graph$as_default(graph)

#load data

setwd("/home/cjls4/feature_vectors/")

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

# which promoters to look at?

# first the H3K27ac_B 

dim(x_train) # 43826 2750 19

dim(train_g4) # 43827 why on earth is there an extra observation here

H3K27ac_B_indices <- which(train_g4$H3K27ac_B!=0)

H3K27ac_B_promoter <- x_train[,,c(H3K27ac_B_indices[-33588])]

preds <- model %>% predict(promoter)

#I dont think that's quite right

# get all the predicted scores, then get the top 100, and find it's equivalent expression score

all_preds <- model %>% predict(x_total)

#check that this is the right column???
top_100_preds_check <- S4Vectors::tail(as.data.table(sort(all_preds[,2] , decreasing = F, index.return = T)), 100)

top_100_preds <- S4Vectors::tail(as.data.table(sort(all_preds[,1] , decreasing = F, index.return = T)), 100)

top_100_names <- x_names[top_100_preds$ix]

top_100_names_check <- x_names[top_100_preds_check$ix]


save(top_100_names, file = "top_100_names.RData")

# bcellf

bf1 <-  total_g4$bcellf[top_100_preds$ix]

# bcellm

bm1 <- total_g4$bcellm[top_100_preds$ix]


# tcellf

tf1 <- total_g4$tcellf[top_100_preds$ix]

# tcellm

tm1 <- total_g4$tcellm[top_100_preds$ix]


# uh corr matrix of b and t cells expression of the top 100 g4s predicted

btmf_df <- cbind(as.integer(bf1), as.integer(bm1), as.integer(tf1), as.integer(tm1))

btmf_df <- as.data.frame(btmf_df)

names(btmf_df) <- c("bf", "bm", "tf", "tm")

save(btmf_df, file = "btmf_df.RData")

heatmap(btmf_df, col = cm.colors(256))
