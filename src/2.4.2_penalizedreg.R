library(MASS)

set.seed(3)

n <- 100
p <- 50

CovMatrix <- outer(1:p, 1:p, function(x,y) {.7^abs(x-y)})

X <- mvrnorm(n, rep(0, p), CovMatrix)
y <- 10*apply(X[,1:2], 1, sum) +
  5*apply(X[,3:5], 1, sum) + 
  apply(X[,6:15], 1, sum) + 
  rnorm(n, 0, 3)
hist(y)
dat = data.frame(cbind(X, y))
names(dat) <- c(paste("X", 1:50, sep = ""), "y")
write.csv(dat, "dat50.csv", row.names = F)

# stepwise
lm.fit <- lm(y~., data = dat)
summary(lm.fit)

null.fit <- lm(y~1, data = dat)
full.fit <- lm(y~., data = dat)

step(null.fit, scope = list(
  lower = null.fit,
  upper = full.fit),
  trace = 0, direction = "forward")
step(full.fit, scope=list(upper=full.fit), trace=0, direction="backward")
step(null.fit, scope = list(upper=full.fit), data=dat, trace=0, direction="both")

## penalized regression

library(glmnet)
train_rows <- sample(1:n, .66*n)
X.train <- X[train_rows,]
X.test <- X[-train_rows,]

y.train <- y[train_rows]
y.test <- y[-train_rows]

fit.lasso <- glmnet(X.train, y.train, family = "gaussian", alpha = 1)
fit.ridge <- glmnet(X.train, y.train, family = "gaussian", alpha = 0)
fit.elnet <- glmnet(X.train, y.train, family = "gaussian", alpha = 0.5)

fit.lasso.cv <- cv.glmnet(X.train, y.train, type.measure="mse", alpha=1, 
                          family="gaussian")
fit.ridge.cv <- cv.glmnet(X.train, y.train, type.measure="mse", alpha=0,
                          family="gaussian")
fit.elnet.cv <- cv.glmnet(X.train, y.train, type.measure="mse", alpha=.5,
                          family="gaussian")

par(mfrow = c(1,3))
plot(fit.lasso, xvar="lambda", main="LASSO")

plot(fit.ridge, xvar="lambda", main="Ridge")

plot(fit.elnet, xvar="lambda",main="Elastic Net")

## logistic regression with penalizing
rm(list = ls())
heart.comp <- read.csv("./Data_Day2and3/heart_comp.csv")

dim(heart.comp)

library(Hmisc)
describe(heart.comp)
heart.fit <- glm(AHD~.,data=heart.comp, family = "binomial")
summary(heart.fit)

library(boot)
cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
LOOCV.err <- cv.glm(heart.comp,heart.fit, cost, K = nrow(heart.comp))$delta
LOOCV.err
cv.10.err <- cv.glm(heart.comp,heart.fit, cost, K = 10)$delta
cv.10.err
library(caret)
data(Default, package = "ISLR")
set.seed(42)
default_idx = createDataPartition(Default$default, p = 0.75, list = FALSE)
default_trn = Default[default_idx, ]
default_tst = Default[-default_idx, ]

cv_5 = trainControl(method = "cv", number = 5)
enet1 = train(
  default ~ ., data = default_trn,
  method = "glmnet",
  trControl = cv_5,
  #tuneLength = 10
)
enet1

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

get_best_result(enet1)

calc_acc = function(actual, predicted) {
  mean(actual == predicted)
}
acc1 <- calc_acc(actual = default_tst$default,
                 predicted = predict(enet1, newdata = default_tst))
print(acc1)

# new data split
default_idx = createDataPartition(Default$default, p = 0.5, list = FALSE)
default_trn = Default[default_idx, ]
default_tst = Default[-default_idx, ]

enet2 = train(
  default ~ ., data = default_trn,
  method = "glmnet",
  trControl = cv_5
  #tuneLength = 10
)

# To simply extract the best result (or tuning parameters) using a quick helper function.
get_best_result(enet2)

# test accuracy using test set
accur2 <- calc_acc(actual = default_tst$default,
                   predicted = predict(enet2, newdata = default_tst))
print(accur2)