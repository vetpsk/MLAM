library(tidyverse)
library(magrittr)
library(MASS)
library(parallel)
library(plotly)
library(caret)
library(glmnet)
install.packages("ranger")
library(ranger)
library(OmicsPLS)

load("DownSyndrome.RData")
str(methylation)
str(glycomics)
str(CpG_groups)
boxplot(glycomics)
boxplot(methylation[,1:100])
print(sample(ClinicalVars, size = 10, replace = T))
table(ClinicalVars$group)

meth <- cbind(methylation, ClinicalVars[,3:5])
set.seed(2021)
inTraining <- createDataPartition(meth$group_ds, p = .70, list = FALSE)
training  <- meth[ inTraining,]
testing  <- meth[-inTraining,]
prop.table(table(training$group_ds))
prop.table(table(testing$group_ds))

cv_5 = trainControl(method = "cv", number = 5)

meth_elnet = train(group_ds ~ ., data = training,
                   method = "glmnet",
                   trControl = cv_5
)

meth_elnet

meth_elnet_int = train(
  group_ds ~ ., data = training,
  method = "glmnet",
  trControl = cv_5,
  tuneLength = 10
)

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

get_best_result(meth_elnet_int)

myControl <- trainControl(
  method = "repeatedcv", number = 10, repeats = 10,
  summaryFunction = twoClassSummary,
  classProbs = TRUE, 
  verboseIter = FALSE,
  savePredictions = TRUE,
)

set.seed(202106)

model1 <- train(group_ds ~ ., training, method = "glmnet", trControl = myControl)
model1


model2 <- train(
  group_ds ~ ., data = training,
  metric = "ROC",
  method = "glmnet",
  tuneGrid = expand.grid(
    alpha = (5:10)/10,
    lambda = 0:10/10
  ),
  trControl = myControl
)
model2
plot(model2)
plot(model2$finalModel)

set.seed(202106)

model_rf <- train(
  group_ds ~ ., training,
  metric = "ROC",
  method = "ranger",
  trControl = myControl
)

model_rf

modelList <- list(
  glmnet1=model1,
  glmnet2 = model2,
  rf = model_rf
)

bwplot(resamples(modelList),
       metric = "ROC")

# Validate on test set with ensemble
allPreds <- sapply(modelList, predict, newdata = testing)

allPreds

# For visualizing the results of grid search, and computing predicted values and cofusion matrix.
ggplot(model1)

model1Class = predict (model1, newdata = testing)
model1Probs <- predict(model1, newdata = testing, type = "prob")



confusionMatrix( as.factor(allPreds[,1]), as.factor(testing$group_ds))
confusionMatrix( as.factor(allPreds[,2]), as.factor(testing$group_ds))
confusionMatrix( as.factor(allPreds[,3]), as.factor(testing$group_ds))

set.seed(83654)
crossval_o2m_adjR2(methylation, glycomics, 1:5, 0:10, 0:9, 10, 8)

par(mfrow=c(1,3))
plot(svd(crossprod(methylation,glycomics),0,0)$d^2 %>% (function(e) e/sum(e)), main='Joint Scree plot')
plot(svd(tcrossprod(methylation),0,0)$d %>% (function(e) e/sum(e)), main="Methylation Scree plot")
plot(svd(crossprod(glycomics),0,0)$d %>% (function(e) e/sum(e)), main="Glycomics Scree plot")
par(mfrow=c(1,1))

r <- 2; rx <- 4; ry <- 6

length(CpG_groups)
length(unique(CpG_groups))

set.seed(7545)
crossval_sparsity(methylation, glycomics, r,rx,ry, 10, 
                  groupx = CpG_groups, keepx_seq = 1:10*40, keepy_seq = 10)


cv_folds <- 10 # number of folds
cv_intervals <- 10 # a grid with keepx values (groups to keep)
cv_set <- c(DS=cut(1:29, cv_folds) %>% as.numeric, # define folds in each outcome class
            SB=cut(1:27, cv_folds) %>% as.numeric, 
            MA=cut(1:29, cv_folds) %>% as.numeric)
X <- methylation
Y <- glycomics
group <- ClinicalVars$group_ds
cv_grid <- round(seq(1,length(unique(CpG_groups)),length.out=cv_intervals))

cvoutp <- mclapply(mc.cores=1,cv_grid, # compute errors across the grid, over the folds
                   function(jj){
                     jj <- as.numeric(jj)
                     cat(jj, "/", length(unique(CpG_groups)), "; ",sep="")
                     outp <- sapply(1:cv_folds, function(ii){
                       iii <- which(ii == cv_set)
                       Xtst <- X[iii,]
                       Ytst <- Y[iii,]
                       Xtrn <- X[-iii,]
                       Ytrn <- Y[-iii,]
                       fit <- try(suppressMessages(o2m(Xtrn,Ytrn,r,rx,ry,sparse=T,keepx=jj,groupx = CpG_groups)), silent = TRUE)
                       if(inherits(fit, "try-error")) fit <- suppressMessages(o2m(Xtrn,Ytrn,r,rx,ry,sparse=T,keepx=jj-1,groupx = CpG_groups))
                       err1 <- mse(Y,predict(fit, X))/sqrt(ssq(Y))
                       err2 <- mse(X,predict(fit, Y, "Y"))/sqrt(ssq(X))
                       fit_outc <- lda(x=fit$Tt,group[-iii])
                       err3 <- mse(as.numeric(group[iii])-1, 
                                   as.numeric(predict(
                                     fit_outc,(Xtst-Xtst%*%fit$W_Y%*%t(fit$P_Y))%*%fit$W.)$class
                                   )-1)/sqrt(ssq(as.numeric(group[iii])-1))
                       c(nr = jj, Yhat = err1, Xhat = err2, outc = err3)
                     })
                     outp
                   })
names(cvoutp) <- as.character(cv_grid)

sapply(cvoutp, function(e) rowMeans(e)) %>% t %>% 
  as.data.frame %>% mutate(across(-nr, scale)) %>% 
  gather(key = "Type", value = "Error", -nr) %>% 
  ggplot(aes(x=nr, y=Error, col=Type)) + geom_line()


cv_gridbest <- round(c(25, 50, 75, 100, 150, 175, 200))
cvoutp_best <- mclapply(mc.cores = 1, cv_gridbest, 
                        function(jj){
                          jj <- as.numeric(jj)
                          cat(jj,";  ")
                          outp <- sapply(1:cv_folds, function(ii){
                            iii <- which(ii == cv_set)
                            Xtst <- X[iii,]
                            Ytst <- Y[iii,]
                            Xtrn <- X[-iii,]
                            Ytrn <- Y[-iii,]
                            fit <- try(suppressMessages(o2m(Xtrn,Ytrn,r,rx,ry,sparse=T,keepx=jj,groupx = CpG_groups)), silent = TRUE)
                            if(inherits(fit, "try-error")) fit <- suppressMessages(o2m(Xtrn,Ytrn,r,rx,ry,sparse=T,keepx=jj-1,groupx = CpG_groups))
                            err1 <- mse(Y,predict(fit, X))/sqrt(ssq(Y))
                            err2 <- mse(X,predict(fit, Y, "Y"))/sqrt(ssq(X))
                            fit_outc <- lda(x=fit$Tt,group[-iii])
                            err3 <- mse(as.numeric(group[iii])-1, 
                                        as.numeric(predict(
                                          fit_outc,(Xtst-Xtst%*%fit$W_Y%*%t(fit$P_Y))%*%fit$W.)$class
                                        )-1)/sqrt(ssq(as.numeric(group[iii])-1))
                            c(nr = jj, Yhat = err1, Xhat = err2, outc = err3)
                          })
                          outp
                        }); 
names(cvoutp_best) <- as.character(cv_gridbest)

lapply(cvoutp_best, function(e) t(e) %>%
         as.data.frame %>% gather(key = "Type", value = "Error",-nr)) %>%
  Reduce(f=bind_rows) %>%
  ggplot(aes(x=as.factor(nr), y=Error)) + geom_boxplot() +
  facet_grid(Type~.,scales = "free")

fit <- o2m(methylation, glycomics, r, rx, ry, sparse = TRUE, 
           groupx = CpG_groups, keepx = 175)
summary(fit)

plot(fit,loading_name = "Yj",i=1,j=2,label = "col")
plot(fit,loading_name = "gr_Xj",i=1,j=2,label = "col", alpha = (apply(fit$W_gr,1,prod) > 0))

plot(data.frame(Tt=fit$Tt,U=fit$U), col = ClinicalVars$group_ds)
data.frame(Group = ClinicalVars$group, JPC = scores(fit, "Xjoint")) %>% 
  pivot_longer(-Group,
               names_to = "Comp", 
               values_to = "LoadingValue") %>% 
  ggplot(aes(x=Comp, y=LoadingValue, col=Group)) + 
  geom_boxplot() + 
  theme_bw()

glm(ClinicalVars$group_ds ~ ., data = data.frame(JPC=scores(fit, "Xjoint")), 
    family = "binomial") %>% 
  summary()



gene_list <- row.names(loadings(fit, loading_name = "gr_Xj", subset = 2))
gene_list %<>% paste0(collapse = ";") %>% 
  str_split(";") %>% unlist %>% unique 
gene_list %>% cat(sep="\n")
