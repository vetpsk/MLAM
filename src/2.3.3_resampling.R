library(MASS)

corX <- matrix(0.5, 8, 8); diag(corX) = 1;
nsize = 100; mu = c(rep(0,8)); covX = corX;
X = mvrnorm(nsize, mu, covX)
mean(X); cov(X)

Betas = c(3, 1.5, 0, 0 ,2, 0, 0, 0)
error <- rnorm(100, 0, 3)
y = X %*% Betas + error
hist(y)

dat <- data.frame(cbind(X, y))
names(dat) <- c(paste("X", 1:8, sep = ""), "y")
write.csv(dat, "dat1.csv", row.names = F)

dat <- read.csv("dat1.csv")

lm.fit <- lm(y ~ ., data = dat)
summary(lm.fit)

err <- y-predict(lm.fit)
plot(err)
qqnorm(err, pch = 1, frame = FALSE)
qqline(err, col = "steelblue", lwd = 2)

heatmap(X)
cor(X)


rm(list = ls())
corX <- matrix(0.5, 8, 8); diag(corX) = 1;
nsize = 10000; mu = c(rep(0,8)); covX = corX;
X = mvrnorm(nsize, mu, covX)
mean(X); cov(X)

Betas = c(3, 1.5, 0, 0 ,2, 0, 0, 0)
error <- rnorm(100, 0, 3)
y = X %*% Betas + error
hist(y)

dat2 <- data.frame(cbind(X, y))
names(dat2) <- c(paste("X", 1:8, sep = ""), "y")
write.csv(dat2, "dat2.csv", row.names = F)

dat <- read.csv("dat2.csv")

lm.fit <- lm(y ~ ., data = dat2)
summary(lm.fit)

err <- y-predict(lm.fit)
plot(err)
qqnorm(err, pch = 1, frame = FALSE)
qqline(err, col = "steelblue", lwd = 2)

heatmap(X)
cor(X)



rm(list = ls())
nsize = 10000
train = sample(1:nsize, 0.66*nsize)
dat <- read.csv("dat2.csv")
lm.fit.1 <- lm(y~.,data =dat, subset = train)
summary(lm.fit.1)
mean((dat$y-predict(lm.fit.1, dat))[-train]^2)
mean((dat$y-predict(lm.fit.1, dat))[train]^2)

lapply(1:10, function(x) {
  train = sample(1:nsize, 0.66*nsize)
  dat <- read.csv("dat2.csv")
  lm.fit.1 <- lm(y~.,data =dat, subset = train)
  mean((dat$y-predict(lm.fit.1, dat))[-train]^2)
})


library(boot)
dat <- read.csv("dat1.csv")
glm.fit = glm(y~., data = dat)
coef(glm.fit)
lm.fit <- lm(y~., data= dat)
coef(lm.fit)
cv.err = cv.glm(dat, glm.fit)
cv.err$delta

set.seed(20200623)
cv.error.5 <- rep(0,5)
for (i in 1:5) {
  glm.fit <- glm(y~., data=dat)
  cv.error.5[i]=cv.glm(dat,glm.fit,K=10)$delta[1]
}
cv.error.5

boot.fn <- function(data, index) return(coef(lm(y~., data=dat, subset = index)))

boot.fn(dat, 1:100)
boot(dat, boot.fn, 1000)
summary(lm(y~.,data=dat))$coef
