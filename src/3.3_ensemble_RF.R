data1 <- readRDS("./Data_Day2and3/car.rds")
head(data1)

str(data1)
summary(data1)

set.seed(100)
train <- sample(nrow(data1), 0.7*nrow(data1), replace = FALSE)
TrainSet <- data1[train,]
ValidSet <- data1[-train,]

library(randomForest)

model1 <- randomForest(Condition ~ ., data = TrainSet, importance = TRUE)
model1
model2 <- randomForest(Condition ~ ., data = TrainSet, ntree = 500, mtry = 6, importance = TRUE)
model2

predTrain <- predict(model2, TrainSet, type = "class")
table(predTrain, TrainSet$Condition)  
predValid <- predict(model2, ValidSet, type = "class")
mean(predValid == ValidSet$Condition) 
table(predValid,ValidSet$Condition)


importance(model2) 
varImpPlot(model2)   
