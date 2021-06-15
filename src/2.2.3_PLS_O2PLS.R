?OmicsPLS
?o2m
library(OmicsPLS)

load("rna_metab.Rdata")

fit <- o2m(rna, metab, 2, 1, 10) #(X, Y, n, nx, ny)
summary(fit)


library(magrittr)
library(ggplot2)
library(gridExtra)
library(illuminaHumanv3.db)
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

######### Plot loadings with OmicsPLS plot method ###
p_metab <- plot(fit, loading_name="Yj", i=1, j=2, label="c", # Plot the loadings
                alpha=0) + # set points to be 100% transparant
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

alprna <- loadings(fit, "Xjoint", 1:2) %>% raise_to_power(2) %>% rowSums
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
p_rna <- ggplot(data.frame(x = fit$W.[, 1], y = fit$W.[, 2]), 
                aes(x = x, y = y),
                alpha = alprna,
                aes(label = NA)) +
  ##################### Add all layers ###
  theme_bw() +
  coord_fixed(.8, c(-.15,.15),c(-.15,.15)) +
  geom_point(alpha = 0.5, col = 'grey') +
  geom_point(data = data.frame(x=fit$W.[LLnr,1],y=fit$W.[LLnr,2]),
             shape = 2, col = 2, size = 2) + 
  geom_text(data = data.frame(x=fit$W.[toprna,1],y=fit$W.[toprna,2]),
            hjust = rep(c(1, 0), length.out = length(toprna)),
            aes(label = names_rna)) + 
  labs(title = "Transcript joint loadings",
       x = "First Joint Loadings", y = "Second Joint Loadings") +
  theme(plot.title = element_text(face='bold')) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

## Finally plot both plots in one figure.
grid.arrange(p_metab, p_rna, ncol=2)


## exercise 2
svd1 <- svd(crossprod(metab, rna), 2, 2)
fit$C.
svd(crossprod(metab, rna), 2, 2)$u
ssq(fit$C. - svd1$u)
