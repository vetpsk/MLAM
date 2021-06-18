##Challenge
##Your task, or challenge, is to build the best possible model that predicts who has heart disease. The hd (0/1) target variable in the training data contains this information. 
##Step 1. You start by reading the two datafiles:

traindata <- readRDS("heart1case.RDS")
testdata <- readRDS("heart2test.RDS")

#Loading packages

library(tree)
library(Hmisc)
library(randomForest)
library(missForest)



#Summarize data
Hmisc::describe(traindata)


###SIMPLE TREE####
# The unpruned tree
tree.Heart=tree(hd~.,traindata) 
summary(tree.Heart)

#The summary() function lists the variables that are used as internal nodes
#in the tree, the number of terminal nodes, and the (training) error rate.
#Q1. What is the training error rate? Note that for classification trees, the deviance
#is reported. A small deviance indicates a tree that provides
#a good fit to the (training) data. 

plot(tree.Heart)
text(tree.Heart,pretty=0, cex=1)

#Getting the test error
set.seed(123)
train=sample(1:nrow(traindata), 170) 
Heart.test=traindata[-train,]
hd.test=Heart.test$hd
tree.Heart=tree(hd~.,traindata,subset=train)
tree.pred=predict(tree.Heart,Heart.test,type="class")
table(tree.pred,hd.test)
#calculate according to the rates interested

##Pruning the tree
#The function cv.tree() performs cross-validation in order to
#determine the optimal level of tree complexity; cost complexity pruning
#is used in order to select a sequence of trees for consideration.
set.seed(123)
cv.Heart=cv.tree(tree.Heart,FUN=prune.misclass)
names(cv.Heart)
cv.Heart
#Note that dev corresponds to the cross-validation error rate 
#in this case. Which is the tree with the lowest cross-validation 
#error rate? Plot the error rate as a function of both size and k
par(mfrow=c(1,2))
plot(cv.Heart$size,cv.Heart$dev,type="b")
plot(cv.Heart$k,cv.Heart$dev,type="b")

set.seed(123)
prune.Heart=prune.misclass(tree.Heart,best=4)
plot(prune.Heart)
text(prune.Heart,pretty=0)
tree.pred=predict(prune.Heart,Heart.test,type="class")
tree.pred
table(tree.pred,hd.test)

###RANDOM FOREST#####
#NRMSE is normalized mean squared error. It is used to 
#represent error derived from imputing continuous values. 
#PFC (proportion of falsely classified) 
#is used to represent error derived from imputing categorical values.

traindata_NA <- traindata[complete.cases(traindata),]
#to much missing values reducing observations from 213 ro 90


##Imputing the data
dat.imp <- missForest(traindata)
dat.imp$OOBerror
traindata.imp <- missForest(traindata)$ximp


model1 <- randomForest(hd ~ ., data = traindata.imp, importance = TRUE)
model1

# Predicting on train set
predTrain <- predict(model1, traindata.imp, type = "class")
# Checking classification accuracy
table(predTrain, traindata.imp$hd)  

# Predicting on Validation set
predValid <- predict(model1, testdata, type = "class")

# Checking classification accuracy
mean(predValid == testdata$hd)                    
table(predValid, testdata$hd)

importance(model1)        
varImpPlot(model1)  












