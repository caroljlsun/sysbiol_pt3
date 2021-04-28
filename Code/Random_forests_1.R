
library(randomForest)
library(tidyverse)
library(caret)
library(e1071)
library(stats)
library(data.table)

setwd("/home/cjls4/feature_vectors/")


#import data
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

#as factor
train_g4$G4 <- as.factor(train_g4$G4)
test_g4$G4 <- as.factor(test_g4$G4)


train_g4.rf <- randomForest(G4~ R_loops + GC, data = train_g4, ntree = 500, proximity = T)
train_g4.rf
plot(train_g4.rf)

# Call:
#   randomForest(formula = G4 ~ R_loops + GC, data = train_g4, ntree = 100,      proximity = T) 
# Type of random forest: classification
# Number of trees: 100
# No. of variables tried at each split: 1
# 
# OOB estimate of  error rate: 31.79%
# Confusion matrix:
#   0    1 class.error
# 0 29858 3567   0.1067165
# 1 10364   37   0.9964426

# uh, some tuning 

model_tuned <- tuneRF(
  x=train_g4[,-20], #define predictor variables
  y=train_g4$G4, #define response variable
  ntreeTry=100,
  mtryStart=4, 
  stepFactor=1.5,
  improve=0.01,
  trace=T, #show real-time progress
  doBest = T,
  plot = T
)

plot(model_tuned)

#Call:
# randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1]) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 2
# 
# OOB estimate of  error rate: 23.09%
# Confusion matrix:
#   0   1 class.error
# 0 32939 486  0.01454001
# 1  9632 769  0.92606480

# names of factors

varNames <- names(train_g4)
varNames <- varNames[!varNames %in% c("G4")]
varNames_1 <- paste(varNames, collapse = "+")
extra.formula <- as.formula(paste("G4", varNames_1, sep = " ~ "))

#methylation_train_g4.rf <- randomForest(G4 ~ Methylation, data = train_g4, ntree = 100, proximity = T)

extra_train_g4.rf <- randomForest(extra.formula,
                                  train_g4,
                                  ntree = 500,
                                  importance = T,
                                  proximity = T)

plot(extra_train_g4.rf)

#extra_train_g4.rf

# Call:
#   randomForest(formula = extra.formula, data = train_g4, ntree = 500,      importance = T, proximity = T) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 4
# 
# OOB estimate of  error rate: 23.65%
# Confusion matrix:
#   0    1 class.error
# 0 31881 1544  0.04619297
# 1  8819 1582  0.84789924

# Plot the importance of different variables

#sort orders your variables from least to most or vs
# order
#main refers to the title of the plot
# n.var = no. of variables to show


varImpPlot(train_g4.rf, 
           sort = T, 
           main = "Variable importance", 
           n.var = 10)

varImpPlot(model_tuned, 
           sort = T, 
           main = "Variable importance", 
           n.var = 10)

varImpPlot(extra_train_g4.rf, 
           sort = T, 
           main = "Variable importance", 
           n.var = 19)


##### predictions #####

### basic model ###
train_g4$predicted.response <- predict(train_g4.rf, train_g4)

#make confusion matrix for the training data

confusionMatrix(data = train_g4$predicted.response, 
                reference = train_g4$G4, 
                positive = "1")
# Accuracy : 0.9992          
# 95% CI : (0.9989, 0.9995)
# No Information Rate : 0.7627          
# P-Value [Acc > NIR] : < 2.2e-16

#same but on the test data

test_g4$predicted.response <- predict(train_g4.rf, 
                                      test_g4)

confusionMatrix(data = test_g4$predicted.response,
                reference = test_g4$G4, 
                positive = "1")
# Accuracy : 0.7032          
# 95% CI : (0.6856, 0.7204)
# No Information Rate : 0.7772          
# P-Value [Acc > NIR] : 1      

### tuned model ###

train_g4$predicted.response.tuned <- predict(model_tuned, train_g4)

#make confusion matrix for the training data

confusionMatrix(data = train_g4$predicted.response.tuned, 
                reference = train_g4$G4, 
                positive = "1")

# Accuracy : 0.8547         
# 95% CI : (0.8513, 0.858)
# No Information Rate : 0.7627         
# P-Value [Acc > NIR] : < 2.2e-16    


#same but on the test data

test_g4$predicted.response.tuned <- predict(model_tuned, 
                                            test_g4)

confusionMatrix(data = test_g4$predicted.response.tuned,
                reference = test_g4$G4, 
                positive = "1")
# Accuracy : 0.7835          
# 95% CI : (0.7675, 0.7989)
# No Information Rate : 0.7772          
# P-Value [Acc > NIR] : 0.2233 


### extra model ###

#same thing but using the extra rf

train_g4$predicted.response.extra <- predict(extra_train_g4.rf, train_g4)

#make confusion matrix for the training data

confusionMatrix(data = train_g4$predicted.response.extra, 
                reference = train_g4$G4, 
                positive = "1")

# Accuracy : 0.9993         
# 95% CI : (0.999, 0.9995)
# No Information Rate : 0.7627         
# P-Value [Acc > NIR] : < 2.2e-16   

#same but on the test data

test_g4$predicted.response.extra <- predict(extra_train_g4.rf, 
                                            test_g4)

confusionMatrix(data = test_g4$predicted.response.extra,
                reference = test_g4$G4, 
                positive = "1")

# Accuracy : 0.7909          
# 95% CI : (0.7751, 0.8061)
# No Information Rate : 0.7772          
# P-Value [Acc > NIR] : 0.04498  


### plots of aurocs ###

um <- test_g4$predicted.response.extra

um[test_g4$predicted.response.extra == 1] <- 1

um[test_g4$predicted.response.extra ==2] <- 0

roc_extra <- roc(response = as.numeric(as.character(test_g4$G4)),
                 predictor = as.numeric(as.character(um)),
                 auc = T)
auc_extra <- auc(roc_extra)

plot(roc_extra)
plot(auc_extra)

roc_svm_test <- roc(response = test_set$G4, predictor =as.numeric(svmPred))

plot(roc_svm_test)





#save random forests
setwd("/home/cjls4/ML/")

#can't save things for some reason lol

save.image(file = "randomforests_1.RData")

save(train_g4.rf, file = "train_g4_rf.RData")
save(model_tuned, file = "model_tuned.RData")
save(extra_train_g4.rf, file = "extra_train_g4_rf.RData")

save(train_g4, file = "train_g4.RData")

#### I need r version 3.6 to get tree, ask russell if that's possible

options(repos='http://cran.rstudio.org')
have.packages <- installed.packages()
cran.packages <- c('devtools','plotrix','randomForest','tree')
to.install <- setdiff(cran.packages, have.packages[,1])
if(length(to.install)>0) install.packages(to.install)

library(devtools)
if(!('reprtree' %in% installed.packages())){
  install_github('araastat/reprtree')
}
for(p in c(cran.packages, 'reprtree')) eval(substitute(library(pkg), list(pkg=p)))

#example tree
library(randomForest)
library(reprtree)

model <- randomForest(Species ~ ., data=iris, importance=TRUE, ntree=500, mtry = 2, do.trace=100)

reprtree:::plot.getTree(model)

##### top 100 G4 scores #####

load("model_tuned.RData")
# get all the predicted scores, then get the top 100, and find it's equivalent expression score

all_preds_rf <- predict(extra_train_g4.rf, train_g4, type = "prob")
#where 1 is the positive class
top_100_preds_rf <- S4Vectors::tail(as.data.table(sort(all_preds_rf[,2] , decreasing = F, index.return = T)), 100)

load("min_g4_positive_names_train.RData")
load("min_g4_negative_names_train.RData")
load("plus_g4_positive_names_train.RData")
load("plus_g4_negative_names_train.RData")

x_names <- rbind(min_g4_positive_names_train, min_g4_negative_names_train, plus_g4_positive_names_train,
                 plus_g4_negative_names_train)

top_100_names_rf <- x_names[top_100_preds_rf$ix]


save(top_100_names_rf, file = "top_100_names_rf.RData")

save(all_preds_rf, file = "all_preds_rf.RData")
