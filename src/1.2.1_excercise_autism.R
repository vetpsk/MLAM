getwd()

dat <- read.csv("Autism.csv")

dat <- dat[dat$verbal == "verbal",]
head(dat)

Items <- dat[,-(1:3)]
dat1 <- dat[,-(1:2)]

pca.Items <- prcomp(Items, scale = TRUE)

summary(pca.Items)

par(mfrow = c(1,1), mai = c(1,1.5,0.2,1.5))

plot(pca.Items$x[,1], pca.Items$x[,2],
     xlab = "PC1", ylab = "PC2",
     col = as.integer(dat$Genetic_Syndrome),
     pch = as.character(as.integer(dat$Genetic_Syndrome)),
     cex = 1.2, cex.axis = 1.3, cex.lab = 1.3, asp = 1)
mn.pc1.Items <- tapply(pca.Items$x[,1],INDEX = dat$Genetic_Syndrome,FUN=mean)
mn.pc2.Items <- tapply(pca.Items$x[,2],INDEX = dat$Genetic_Syndrome,FUN=mean)
points(x=mn.pc1.Items,y=mn.pc2.Items,col="black",pch=19,cex=3)
points(x=mn.pc1.Items,y=mn.pc2.Items,col=1:8,pch=19,cex=2)

#supervised learning: binary classification

select14 <- dat1$Genetic_Syndrome %in% levels(dat1$Genetic_Syndrome)[c(1,4)]
table(dat[select14,]$Genetic_Syndrome)

dat.pca <- cbind(dat, pca.Items$x)

datglm.fit <- glm(Genetic_Syndrome == "22q11DS" ~
                    PC1 + PC2, 
                  data = dat.pca[select14,],
                  family = "binomial")

summary(datglm.fit)
# cross-val of this model
cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)

library(boot) #this line is not generally needed as the package is preloaded

cv.err = cv.glm(dat.pca[select14,], datglm.fit, cost = cost)

c(out.of.sample.error = cv.err$delta[2])
in.sample.error <- cost(datglm.fit$fitted.values,
                        dat.pca[select14,]$Genetic_Syndrome=="22q11DS")
c(in.sample.error=in.sample.error)

# 1-NN classification
library(class)

NN1.train = knn(train = dat.pca[select14, c("PC1", "PC2")],
                test =dat.pca[select14,c("PC1","PC2")], 
                cl=dat.pca[select14,c("Genetic_Syndrome")], 
                k = 1, l = 0, prob = FALSE, use.all = TRUE)
table(NN1.train, dat.pca[select14,c("Genetic_Syndrome")])
c(in.sample.error.1NN=
    cost(NN1.train!=dat.pca[select14,c("Genetic_Syndrome")]))

#LOOCV
NN1.cv1 <- knn.cv(train=dat.pca[select14,c("PC1","PC2")], 
                  cl=dat.pca[select14,c("Genetic_Syndrome")], 
                  k = 1, l = 0, prob = FALSE, use.all = TRUE)
table(NN1.cv1,dat.pca[select14,c("Genetic_Syndrome")])
c(out.sample.error.1NN=
    cost(NN1.cv1!=dat.pca[select14,c("Genetic_Syndrome")])) 

########################## Q1 #############################

knn2 <- knn(train=dat.pca[select14,c("PC1","PC2")], 
            test =dat.pca[select14,c("PC1","PC2")], 
            cl=dat.pca[select14,c("Genetic_Syndrome")], 
            k = 2, l = 0, prob = FALSE, use.all = TRUE)
table(knn2,dat.pca[select14,c("Genetic_Syndrome")])
knn2.cv1 <- knn.cv(train=dat.pca[select14,c("PC1","PC2")], 
                  cl=dat.pca[select14,c("Genetic_Syndrome")], 
                  k = 2, l = 0, prob = FALSE, use.all = TRUE)
table(knn2.cv1,dat.pca[select14,c("Genetic_Syndrome")])
c(out.sample.error.2NN=cost(knn2.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))


knn3 <- knn(train=dat.pca[select14,c("PC1","PC2")], 
            test =dat.pca[select14,c("PC1","PC2")], 
            cl=dat.pca[select14,c("Genetic_Syndrome")], 
            k = 3, l = 0, prob = FALSE, use.all = TRUE)
table(knn3,dat.pca[select14,c("Genetic_Syndrome")])
knn3.cv1 <- knn.cv(train=dat.pca[select14,c("PC1","PC2")], 
                   cl=dat.pca[select14,c("Genetic_Syndrome")], 
                   k = 3, l = 0, prob = FALSE, use.all = TRUE)
table(knn3.cv1,dat.pca[select14,c("Genetic_Syndrome")])
c(out.sample.error.2NN=cost(knn3.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))


knn10 <- knn(train=dat.pca[select14,c("PC1","PC2")], 
            test =dat.pca[select14,c("PC1","PC2")], 
            cl=dat.pca[select14,c("Genetic_Syndrome")], 
            k = 10, l = 0, prob = FALSE, use.all = TRUE)
table(knn10,dat.pca[select14,c("Genetic_Syndrome")])
knn10.cv1 <- knn.cv(train=dat.pca[select14,c("PC1","PC2")], 
                   cl=dat.pca[select14,c("Genetic_Syndrome")], 
                   k = 10, l = 0, prob = FALSE, use.all = TRUE)
table(knn10.cv1,dat.pca[select14,c("Genetic_Syndrome")])
c(out.sample.error.2NN=cost(knn10.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))

knn15 <- knn(train=dat.pca[select14,c("PC1","PC2")], 
             test =dat.pca[select14,c("PC1","PC2")], 
             cl=dat.pca[select14,c("Genetic_Syndrome")], 
             k = 15, l = 0, prob = FALSE, use.all = TRUE)
table(knn15,dat.pca[select14,c("Genetic_Syndrome")])
knn15.cv1 <- knn.cv(train=dat.pca[select14,c("PC1","PC2")], 
                    cl=dat.pca[select14,c("Genetic_Syndrome")], 
                    k = 15, l = 0, prob = FALSE, use.all = TRUE)
table(knn15.cv1,dat.pca[select14,c("Genetic_Syndrome")])
c(out.sample.error.2NN=cost(knn15.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))

knn7 <- knn(train=dat.pca[select14,c("PC1","PC2")], 
             test =dat.pca[select14,c("PC1","PC2")], 
             cl=dat.pca[select14,c("Genetic_Syndrome")], 
             k = 7, l = 0, prob = FALSE, use.all = TRUE)
table(knn7,dat.pca[select14,c("Genetic_Syndrome")])
knn7.cv1 <- knn.cv(train=dat.pca[select14,c("PC1","PC2")], 
                    cl=dat.pca[select14,c("Genetic_Syndrome")], 
                    k = 7, l = 0, prob = FALSE, use.all = TRUE)
table(knn7.cv1,dat.pca[select14,c("Genetic_Syndrome")])
c(out.sample.error.2NN=cost(knn7.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))

c(out.sample.error.1NN=cost(NN1.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))
c(out.sample.error.2NN=cost(knn2.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))
c(out.sample.error.3NN=cost(knn3.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))
c(out.sample.error.7NN=cost(knn7.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))
c(out.sample.error.10NN=cost(knn10.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))
c(out.sample.error.15NN=cost(knn15.cv1!=dat.pca[select14,c("Genetic_Syndrome")]))

########################## Q2 #############################

select15 <- dat1$Genetic_Syndrome %in% levels(dat1$Genetic_Syndrome)[c(1,5)]
table(dat[select15,]$Genetic_Syndrome)

datglm.fit2 <- glm(Genetic_Syndrome == "22q11DS" ~
                    PC1 + PC2, 
                  data = dat.pca[select15,],
                  family = "binomial")

summary(datglm.fit2)
# cross-val of this model
cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)

library(boot) #this line is not generally needed as the package is preloaded

cv.err2 = cv.glm(dat.pca[select15,], datglm.fit2, cost = cost)

c(out.of.sample.error = cv.err2$delta[2])
in.sample.error <- cost(datglm.fit2$fitted.values,
                        dat.pca[select15,]$Genetic_Syndrome=="22q11DS")
c(in.sample.error=in.sample.error)

# 1-NN classification
library(class)

NN1.train2 = knn(train = dat.pca[select15, c("PC1", "PC2")],
                test =dat.pca[select15,c("PC1","PC2")], 
                cl=dat.pca[select15,c("Genetic_Syndrome")], 
                k = 1, l = 0, prob = FALSE, use.all = TRUE)
table(NN1.train2, dat.pca[select15,c("Genetic_Syndrome")])
c(in.sample.error.1NN2=
    cost(NN1.train2!=dat.pca[select15,c("Genetic_Syndrome")]))

#LOOCV
NN1.cv2 <- knn.cv(train=dat.pca[select15,c("PC1","PC2")], 
                  cl=dat.pca[select15,c("Genetic_Syndrome")], 
                  k = 1, l = 0, prob = FALSE, use.all = TRUE)
table(NN1.cv2,dat.pca[select15,c("Genetic_Syndrome")])
c(out.sample.error.1NN2=
    cost(NN1.cv2!=dat.pca[select15,c("Genetic_Syndrome")])) 

# lapply to loop over multiple knn models, saves time!  
knn_mods <-  lapply(1:10, function(i) {
  knn(train=dat.pca[select15,c("PC1","PC2")], 
             test =dat.pca[select15,c("PC1","PC2")], 
             cl=dat.pca[select15,c("Genetic_Syndrome")], 
             k = i, l = 0, prob = FALSE, use.all = TRUE)
})

lapply(knn_mods, 
       function(i) c(in.sample.error=
                                 cost(i!=dat.pca[select15,c("Genetic_Syndrome")])))  

