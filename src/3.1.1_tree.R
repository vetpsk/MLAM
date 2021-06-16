library(tree)

Heart <- read.csv("./Data_Day2and3/heart_all.csv",
                  sep = ",", header = T, as.is = F)

library(Hmisc)
describe(Heart)

?tree::tree
tree.Heart <- tree(AHD~., Heart)
summary(tree.Heart)
plot(tree.Heart)
text(tree.Heart, pretty = 0, cex = 0.7, pos = 2)

#test error
set.seed(20200624)
train=sample(1:nrow(Heart), 153) 
Heart.test=Heart[-train,]
AHD.test=Heart.test$AHD
tree.Heart=tree(AHD~.,Heart,subset=train)
tree.pred=predict(tree.Heart,Heart.test,type="class")

table(tree.pred,AHD.test)

set.seed(24062020)
cv.Heart=cv.tree(tree.Heart,FUN=prune.misclass)
names(cv.Heart)
cv.Heart

par(mfrow=c(1,2))
plot(cv.Heart$size,cv.Heart$dev,type="b")
plot(cv.Heart$k,cv.Heart$dev,type="b")

set.seed(17)
prune.Heart=prune.misclass(tree.Heart,best=9)
plot(prune.Heart)
text(prune.Heart,pretty=0)
tree.pred=predict(prune.Heart,Heart.test,type="class")
tree.pred
table(tree.pred,AHD.test)

prune.Heart.12=prune.misclass(tree.Heart,best=12)
plot(prune.Heart.12)
text(prune.Heart.12,pretty=0)
tree.pred=predict(prune.Heart.12,Heart.test,type="class")
tree.pred
table(tree.pred,AHD.test)

cat(c("TP = 44", "TN = 61", "FP = 24", "FN = 21", "\n *considering yes = positive"))
cat(c(paste("Sens = ", 44/(44+21)), paste("Spec = ", 61/(61+24))))
par(mfrow = c(1,1))

rm(list = ls())

# Boston data
library(MASS)

set.seed(1)
train = sample(1:nrow(Boston), nrow(Boston)/2)
tree.boston=tree(medv~.,Boston,subset=train)
summary(tree.boston)
plot(tree.boston)
text(tree.boston, pretty = 0)

cv.boston=cv.tree(tree.boston)
plot(cv.boston$size,cv.boston$dev,type='b')
prune.boston=prune.tree(tree.boston,best=5)
plot(prune.boston)
text(prune.boston,pretty=0)
yhat=predict(tree.boston,newdata=Boston[-train,])
boston.test=Boston[-train,"medv"]
plot(yhat,boston.test)
abline(0,1)
mean((yhat-boston.test)^2)