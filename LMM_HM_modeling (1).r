#Author: trust Odia
#Date Nov 26, 2020
#Objective: implement LMM on repeated measures data

#Run on a HPC
setwd("/home/todia/R_modeling")
load("/home/todia/R_modeling/PETModelling.RData")
MSLAxGenexCell = readRDS("~/R_modeling/MSLAxGenexCell_A.rds")
library(lme4)
library(tidyverse)
library(car)
lmm1.out =list(16147)
lmm1 = list()
lmm2.out =list(16147)
petvsgene = as.data.frame.matrix(matrix(nrow = 16147, ncol = 2))
j = 1
for(i in 5:16151){
  b <- colnames(MSLAxGenexCell)[i]
  lmm1[[b]] = lmer(formula = formula(paste("MSLA.Value ~ ", paste(b, "+ Time + CD8_positive_alpha_beta_T_cell + CD14_positive_monocyte + neutrophil + (1 + Time | SUBJID)"))), data = MSLAxGenexCell_A)
  lmm1.out[[b]] = anova(lmm1[[b]], test = "F")
  lmm1.out[[b]] = lmm1.out[[b]][1,] %>% unlist()
  petvsgene[j,1] = b
  j = j + 1
}
#unpack the list (lmm1.out) into a vector
vec1 = unlist(lmm1.out)
#coerce vector into a dataframe
petvsgenes.1 = tibble(Name = names(vec1), Value = vec1)
#Separate "Name" column
petvsgenes.1 = separate(data = petvsgenes.1, col = Name, into = c("GeneID", "Stat"))
#Spread dataframe
petvsgenes.1 = petvsgenes.1 %>% group_by(Stat) %>% mutate(ind = row_number()) %>% spread(Stat, Value) %>% as.data.frame()
petvsgenes.1 = (petvsgenes.1)[-1,]
petvsgenes.1$`<NA>` = NULL
saveRDS(petvsgenes.1, file = "~/R_modeling/PET.Gene.CellTypes.NoBcell.rds")

