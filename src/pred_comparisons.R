# pranav.pred <- read.csv("predictions_testdata_penreg.csv")
# panos.pred <- read.csv("./src/Test_pred_panos.csv")
# RF.pred <- readxl::read_excel("src/RF.xlsx", col_names = F)
# RF.pred$rf.hd <- RF.pred$...2
# RF.pred$...2 <- NULL
# merge.pred <- cbind(pranav.pred, panos.pred, RF.pred)
# merge.pred[,8] <- NULL
# merge.pred1 <- merge.pred %>% dplyr::select(X, elnet.pred = enet.hd, svm.pred = predict.SVM_model_final..testdata., rf.pred = rf.hd)




merge.pred <- read.csv("gr7_merge_pred.csv")


table(merge.pred$elnet.pred, merge.pred$svm.pred)

table(merge.pred$elnet.pred == merge.pred$svm.pred)


table(merge.pred$elnet.pred, merge.pred$rf.pred)

table(merge.pred$elnet.pred == merge.pred$rf.pred)


table(merge.pred$svm.pred, merge.pred$rf.pred)

table(merge.pred$svm.pred == merge.pred$rf.pred)