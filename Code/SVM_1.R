
library(e1071)
library(caTools)
library(caret)
library(stats)
library(useful)
library(plyr)
library(doMC)
library(pROC)
library(data.table)

# global options

registerDoMC(detectCores()/2)
getDoParWorkers()

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


##### radial kernel function #####

svmRadialE1071 <- list(
  label = "Support Vector Machines with Radial Kernel - e1071",
  library = "e1071",
  type = c("Regression", "Classification"),
  parameters = data.frame(parameter="cost",
                          class="numeric",
                          label="Cost"),
  grid = function (x, y, len = NULL, search = "grid") 
  {
    if (search == "grid") {
      out <- expand.grid(cost = 2^((1:len) - 3))
    }
    else {
      out <- data.frame(cost = 2^runif(len, min = -5, max = 10))
    }
    out
  },
  loop=NULL,
  fit=function (x, y, wts, param, lev, last, classProbs, ...) 
  {
    if (any(names(list(...)) == "probability") | is.numeric(y)) {
      out <- e1071::svm(x = as.matrix(x), y = y, kernel = "radial", 
                        cost = param$cost, ...)
    }
    else {
      out <- e1071::svm(x = as.matrix(x), y = y, kernel = "radial", 
                        cost = param$cost, probability = classProbs, ...)
    }
    out
  },
  predict = function (modelFit, newdata, submodels = NULL) 
  {
    predict(modelFit, newdata)
  },
  prob = function (modelFit, newdata, submodels = NULL) 
  {
    out <- predict(modelFit, newdata, probability = TRUE)
    attr(out, "probabilities")
  },
  predictors = function (x, ...) 
  {
    out <- if (!is.null(x$terms)) 
      predictors.terms(x$terms)
    else x$xNames
    if (is.null(out)) 
      out <- names(attr(x, "scaling")$x.scale$`scaled:center`)
    if (is.null(out)) 
      out <- NA
    out
  },
  tags = c("Kernel Methods", "Support Vector Machines", "Regression", "Classifier", "Robust Methods"),
  levels = function(x) x$levels,
  sort = function(x)
  {
    x[order(x$cost), ]
  }
)


#seperate G4 from the predictors

seg_train_g4 <- train_g4[,1:19]

seg_test_g4 <- test_g4[,1:19]


#preprocessing

transformations <- preProcess(train_g4, 
                              method=c("YeoJohnson", "center", "scale", "corr"),
                              cutoff=0.75)

training_set <- predict(transformations, train_g4)

table(training_set[,20])
table(train_g4[,20])

class(training_set[,20])

training_set$G4 <- as.numeric(training_set$G4)
training_set$G4 <- as.factor(training_set$G4)
training_set$G4 <- as.numeric(training_set$G4)
training_set$G4 <- as.factor(training_set$G4) #yes this is convoluted, but it gets the job done
training_set$G4 <- revalue(training_set$G4, c("1"="A", "2"= "B"))


test_set <- predict(transformations, test_g4)

test_set$G4 <- as.numeric(test_set$G4)
test_set$G4 <- as.factor(test_set$G4)
test_set$G4 <- as.numeric(test_set$G4)
test_set$G4 <- as.factor(test_set$G4) #yes this is convoluted, but it gets the job done
test_set$G4 <- revalue(test_set$G4, c("1"="A", "2"= "B"))



#model cross-validation and tuning

set.seed(42)

seeds <- vector(mode = "list", length = 26)
for(i in 1:25) seeds[[i]] <- sample.int(1000, 9)
seeds[[26]] <- sample.int(1000,1)

# setting cross validation method, trying to tune cost

cvCtrl_probs <- trainControl(method = "repeatedcv", 
                       repeats = 5,
                       number = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       seeds=seeds)

cvCtrl <- trainControl(method = "repeatedcv", 
                       repeats = 5,
                       number = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       seeds=seeds)

# training data
svmTune <- train(x = training_set[,1:19],
                 y = training_set$G4,
                 method = svmRadialE1071,
                 tuneLength = 9,
                 metric = "ROC",
                 trControl = cvCtrl)

probs_svm <- train(x = training_set[,1:19],
                   y = training_set$G4,
                   method = svmRadialE1071,
                   tuneLength = 1,
                   metric = "ROC",
                   trControl = cvCtrl_probs,
                   probability = TRUE)


save(svmTune, file = "svmTune.RData")
save(probs_svm, file = "probs_svm.RData")

svmTune

svmTune$finalModel

plot(svmTune, metric = "ROC", scales = list(x = list(log =2)))

svmPred <- predict(svmTune, test_set[,1:19])

confusionMatrix(svmPred, as.factor(test_set$G4))

# Accuracy : 0.785          
# 95% CI : (0.769, 0.8003)
# No Information Rate : 0.7772         
# P-Value [Acc > NIR] : 0.1717 


##### plot fun things #####

gridSize <- 150 

v1limits <- c(min(test_set$R_loops),max(test_set$R_loops))
tmpV1 <- seq(v1limits[1],v1limits[2],len=gridSize)


v2limits <- c(min(test_set$GC), max(test_set$GC))
tmpV2 <- seq(v2limits[1],v2limits[2],len=gridSize)

xgrid <- expand.grid(tmpV1,tmpV2)
names(xgrid) <- names(training_set)[c(13, 19)]

V3 <- as.numeric(predict(svmTune, xgrid))

V3 <- predict(svmTune, xgrid)


xgrid <- cbind(xgrid, V3)


point_shapes <- c(15,17)
point_colours <- brewer.pal(3,"Dark2")
point_size = 2

trainClassNumeric <- ifelse(moonsTrain$V3=="A", 1, 2)
testClassNumeric <- ifelse(moonsTest$V3=="A", 1, 2)

ggplot(xgrid, aes(V1,V2)) +
  geom_point(col=point_colours[V3], shape=16, size=0.3) +
  geom_point(data=moonsTrain, aes(V1,V2), col=point_colours[trainClassNumeric],
             shape=point_shapes[trainClassNumeric], size=point_size) +
  geom_contour(data=xgrid, aes(x=V1, y=V2, z=V3), breaks=1.5, col="grey30") +
  ggtitle("train") +
  theme_bw() +
  theme(plot.title = element_text(size=25, face="bold"), axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))

ggplot(xgrid, aes(V1,V2)) +
  geom_point(col=point_colours[V3], shape=16, size=0.3) +
  geom_point(data=moonsTest, aes(V1,V2), col=point_colours[testClassNumeric],
             shape=point_shapes[testClassNumeric], size=point_size) +
  geom_contour(data=xgrid, aes(x=V1, y=V2, z=V3), breaks=1.5, col="grey30") +
  ggtitle("test") +
  theme_bw() +
  theme(plot.title = element_text(size=25, face="bold"), axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))

# might need to do dimension reduction etc. to visualise what is going on lol

#oh well. Here's another SVM


training_set[,1:4] = scale(training_set[,1:4])
test_set[,1:4] = scale(test_set[,1:4])


classifier1 = svm(formula = G4~., data = training_set, type = 'C-classification', kernel = 'radial')
classifier2 = svm(formula = G4~ R_loops + GC, data = training_set, type = 'C-classification', kernel = 'radial')

test_pred1 = predict(classifier1, type = 'response', newdata = test_set[,-20])
test_pred2 = predict(classifier2, type = 'response', newdata = test_set[,-20])

# Making Confusion Matrix
cm1 = base::table(unlist(test_set[,20]), test_pred1)
cm2 = table(unlist(test_set[,20]), test_pred2)
cm1 # Confusion Matrix for all parameters
cm2 # Confusion Matrix for parameters being R loops and GsC content



svmPred <- predict(svmTune, test_set[,1:19])

confusionMatrix(test_pred1, as.factor(test_set$G4))
confusionMatrix(test_pred2, as.factor(test_set$G4))
confusionMatrix(test_pred3, as.factor(test_set$G4))

# The accuracy for both model looks solid...

m2 <- svm(Species~., data = iris)

plot(m2, iris, Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))

plot(classifier1, training_set, R_loops ~ GC)

plot(classifier2, training_set, R_loops ~ GC)

iris.part = subset(iris, Species != 'setosa')
iris.part$Species = factor(iris.part$Species)
#iris.part = iris.part[, c(1,2,5)]
svm.fit = svm(formula=Species~., data=iris.part, type='C-classification', kernel='linear')
plot(svm.fit, iris.part, Petal.Width ~ Petal.Length, slice = list(Sepal.Width = 3, Sepal.Length = 4))



##### top 100 G4 scores #####


# get all the predicted scores, then get the top 100, and find it's equivalent expression score

load("svmTune.RData")
load("probs_svm.RData")


varImp(object = svmTune)

roc_svm_test <- roc(response = test_set$G4, predictor =as.numeric(svmPred))

plot(roc_svm_test)




plot(varImp(svmTune), col = viridis(19, direction = 1))
#simple_svm <- svm(x = training_set[,1:19],
                  # y = training_set$G4,
                  # kernel = "radial",
                  # cost = 0.25,
                  # cross = 5,
                  # probability = TRUE)

svmTune$finalModel

#pretty sure "B" is G4 positive

all_preds_svm <-  predict.train(svmTune, type = "prob")

svmPred <- predict(svmTune, test_set[,1:19], type = "prob")

total_preds <- rbind(all_preds_svm, svmPred)

top_100_preds_svm <- S4Vectors::tail(as.data.table(sort(total_preds$B , decreasing = F, index.return = T)), 100)

load("min_g4_positive_names_train.RData")
load("min_g4_negative_names_train.RData")
load("plus_g4_positive_names_train.RData")
load("plus_g4_negative_names_train.RData")
load("min_g4_positive_names_test.RData")
load("min_g4_negative_names_test.RData")
load("plus_g4_positive_names_test.RData")
load("plus_g4_negative_names_test.RData")

x_names <- rbind(min_g4_positive_names_train, min_g4_negative_names_train,
                 plus_g4_positive_names_train, plus_g4_negative_names_train,
                 min_g4_positive_names_test, min_g4_negative_names_test,
                 plus_g4_positive_names_test, plus_g4_negative_names_test)

top_100_names_svm <- x_names[top_100_preds_svm$ix]


save(top_100_names_svm, file = "top_100_names_svm.RData")
