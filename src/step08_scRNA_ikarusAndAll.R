library(pROC)
library(MuSiC)
library(Biobase)
library(Seurat)

source('fwp.R')
source('fwp_source.R')

TYPE1=c('BRCA','KIRP','LIHC','PCPG',
        'BLCA','CESC','CHOL','COAD',
        'ESCA','HNSC','KICH','LUSC',
        'PAAD','THYM','UCEC')

TYPE2=c('scBRCA')
DIR='./data/scRNA/'


OUT=.testAll_scRNA(TYPE1, TYPE2, DIR, .test_ikarus)
saveRDS(OUT, paste0(DIR,'/','eval_ikarus.rds'))




DIR='./data/scRNA/'
O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwo.rds'))
O3=readRDS(paste0(DIR,'/','eval_music.rds'))
O4=readRDS(paste0(DIR,'/','eval_ikarus.rds'))

OOO=cbind(O1[,4],O2[,4],O3[,4],O4[,4])
rownames(OOO)=O1[,1]

pdf('plot/p07_boxplot_scRNA.pdf',width=3,height=3.5)
COL=c(rep('indianred2',1),'royalblue1','grey70','grey70')
boxplot(OOO, pch='+',col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
barplot(OOO['BRCA',],col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
barplot(OOO['CESC',],col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
dev.off()


apply(OOO,2,mean)
#[1] 0.8268625 0.7980829 0.6744619 0.7307102

print(OOO[1,])
#0.9584225 0.9465386 0.8019542 0.5026977


print(OOO['CESC',])
#0.9221674 0.8895294 0.5770416 0.7680788

#plot(O1[,4],O2[,4])
#abline(a=0,b=1)

#t.test(O1[,4],O2[,4],paired=T,alternative='greater')







