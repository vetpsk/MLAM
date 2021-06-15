#2.1.1 PCA exercise

x1 <- rnorm(100)
x2 <- rnorm(100)

X <- cbind(x1, x2)

svd(X, 0, 1)$v

x2 <- x1 + rnorm(100, sd = 0.1)

X2 <- cbind(x1,x2)

svd(X2, 0, 1)$v


## DILGOM study RNA metab data exercise PCA
rm(list = ls())
load("rna_metab.Rdata")
ls()

library(gplots)

heatmap.2(cor(metab, use = "pair"),
          dendrogram = "none",
          Rowv = F, Colv = F, trace = "n",
          breaks = seq(-1, 1, length.out = 25),
          col = gplots::bluered)

par(mfrow = c(2, 1))
boxplot(rna[,1:100])
boxplot(metab)
par(mfrow = c(1,1))

pca.rna <- svd(rna, 0, 2)
pca.metab <- svd(metab, 0, 2)
par(mfrow=c(1,2))
plot(rna %*% pca.rna$v, 
     main = "RNA PCA plot of the scores",
     xlab=NA,ylab=NA)
plot(metab %*% pca.metab$v, 
     main = "Metabolites PCA plot of the scores",
     xlab=NA,ylab=NA)
par(mfrow = c(1,1))

### plot weights

library(magrittr)
library(ggplot2)
library(gridExtra)
library(OmicsPLS)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("illuminaHumanv3.db")
library(illuminaHumanv3.db) # BiocManager::install("illuminaHumanv3.db")
# Color names
LLmodule <- c("ILMN_1690209",'ILMN_1766551', 'ILMN_1749131', 'ILMN_1688423', 
              'ILMN_2102670', 'ILMN_1792323', 'ILMN_1899034', 'ILMN_1806721', 
              'ILMN_1695530', 'ILMN_1726114', 'ILMN_1751625', 'ILMN_1726114', 
              'ILMN_1753648', 'ILMN_1779043')
LLnr <- which(colnames(rna) %in% LLmodule)
rna_genenames <- select(illuminaHumanv3.db, 
                        keys = colnames(rna)[LLnr], 
                        keytype = "PROBEID", columns = "SYMBOL")[,2]

name_col <- 1 + sapply( #First sapply loops over column names
  X = colnames(metab),
  FUN = function(arg){
    crossprod(
      c(1, 1, 3, 4, 5), # Weights to be used as categories
      sapply(c("VLDL", "LDL", "IDL", "HDL","FA"), # metabolite classes
             function(arg2){grepl(arg2, arg)} # compare class of metabolites
      )
    )
  }
)
name_col <- factor(name_col, 
                   levels = c(3,2,4:6,1), 
                   labels = c("VLDL", "LDL", "IDL", "HDL","FA","Other"))

# alpmetab <- loadings(fit, "Yjoint", 1:2) %>%  # Retreive loadings
#   abs %>% # Absolute loading values for positive weights
#   rowSums %>% # Sum over the components
#   sqrt + (name_col!="Other") # Take square root

######### Plot loadings with ggplot ###
p_metab <- ggplot(data.frame(x = pca.metab$v[,1], y = pca.metab$v[, 2]), aes(x = x, y = y)) + 
  ##################### Add all layers ###
  theme_bw() +
  coord_fixed(ratio = 1, xlim=c(-.2,.2),ylim=c(-.2,.2)) +
  geom_point( # Set color and size
    aes(col=name_col, size = I(1+(name_col%in%c("VLDL","HDL"))), 
        shape = name_col),show.legend = T) +
  theme(legend.position="right") +
  scale_color_discrete(name="Metabolite\nGroup",
                       labels=c("VLDL", "LDL", "IDL", "HDL","FA","Other")) +
  guides(size=F) + scale_shape_discrete(name="Metabolite\nGroup",
                                        labels=c("VLDL", "LDL", "IDL", "HDL","FA","Other")) +
  scale_shape_manual(name="Metabolite\nGroup", values=c(15,3,4,17,5,6)) + 
  labs(title = "Metabolite joint loadings",
       x = "First Joint Loadings", y = "Second Joint Loadings") +
  theme(plot.title = element_text(face='bold'),
        legend.title=element_text(face='bold')) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

alprna <- pca.rna$v %>% raise_to_power(2) %>% rowSums
alprna[-(order(alprna,decreasing=T)[1:10])] = 0
alprna <- sign(alprna)
toprna <- which(alprna>0)
names_rna <- mapIds(illuminaHumanv3.db, 
                    keys = colnames(rna)[toprna], 
                    keytype = "PROBEID", 
                    column = "SYMBOL",
                    multiVals = 'first')
names_rna[which(is.na(names_rna))] <- "?"
######### Plot loadings with OmicsPLS plot method ###
p_rna <- ggplot(data.frame(x = pca.rna$v[, 1], y = pca.rna$v[, 2]), 
                aes(x = x, y = y),
                alpha = alprna,
                aes(label = NA)) +
  ##################### Add all layers ###
  theme_bw() +
  coord_fixed(.8, c(-.15,.15),c(-.15,.15)) +
  geom_point(alpha = 0.5, col = 'grey') +
  geom_point(data = data.frame(x = pca.rna$v[LLnr, 1], y = pca.rna$v[LLnr, 2]),
             shape = 2, col = 2, size = 2) + 
  geom_text(data = data.frame(x=pca.rna$v[toprna,1],y=pca.rna$v[toprna,2]),
            hjust = rep(c(1, 0), length.out = length(toprna)),
            aes(label = names_rna)) + 
  labs(title = "Transcript joint loadings",
       x = "First Joint Loadings", y = "Second Joint Loadings") +
  theme(plot.title = element_text(face='bold')) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

## Finally plot both plots in one figure.
grid.arrange(p_metab, p_rna, ncol=2)

### other exercises
plot(pca.rna$v)
identify(pca.rna$v, labels = colnames(rna))
# 2 1793 4103 4819 6121 6211