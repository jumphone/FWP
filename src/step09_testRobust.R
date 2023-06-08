########################
library(pROC)
library(MuSiC)
library(Biobase)
library(Seurat)

source('fwp.R')
source('fwp_source.R')

TYPE=c('BRCA','KIRP','LIHC','PCPG')
DIR='./data/first/'

#OUT=.testAll(TYPE, DIR, .test_fwp, npcs=2)
#saveRDS(OUT, paste0(DIR,'/','eval_fwp.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.9)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.9.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.8)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.8.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.7)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.7.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.6)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.6.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.5)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.5.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.4)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.4.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.3)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.3.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.2)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.2.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.1)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.1.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.01)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.01.rds'))
OUT=.testAll(TYPE, DIR, .test_fwp_dw, npcs=2, dw=0.001)
saveRDS(OUT, paste0(DIR,'/','eval_fwp_dw0.001.rds'))

DIR='./data/first/'
O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwp_dw0.9.rds'))
O3=readRDS(paste0(DIR,'/','eval_fwp_dw0.8.rds'))
O4=readRDS(paste0(DIR,'/','eval_fwp_dw0.7.rds'))
O5=readRDS(paste0(DIR,'/','eval_fwp_dw0.6.rds'))
O6=readRDS(paste0(DIR,'/','eval_fwp_dw0.5.rds'))
O7=readRDS(paste0(DIR,'/','eval_fwp_dw0.4.rds'))
O8=readRDS(paste0(DIR,'/','eval_fwp_dw0.3.rds'))
O9=readRDS(paste0(DIR,'/','eval_fwp_dw0.2.rds'))
O10=readRDS(paste0(DIR,'/','eval_fwp_dw0.1.rds'))
O11=readRDS(paste0(DIR,'/','eval_fwp_dw0.01.rds'))
O12=readRDS(paste0(DIR,'/','eval_fwp_dw0.001.rds'))



OOO=cbind(O1[,4],O2[,4],O3[,4],O4[,4],O5[,4],O6[,4],O7[,4],O8[,4],O9[,4],O10[,4],O11[,4],O12[,4])
apply(OOO,2,mean)

pdf('plot/p10_boxplot_robustDw.pdf',width=6,height=5)
COL=c(rep('indianred2',ncol(OOO)))
boxplot(OOO, pch='+',col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
dev.off()

apply(OOO,2,mean)

#[1] 0.9262461 0.9259464 0.9235663 0.9274571 0.9314842 0.9188348 0.9177012
#[8] 0.9125350 0.9285817 0.9045167 0.8543765 0.6481015




