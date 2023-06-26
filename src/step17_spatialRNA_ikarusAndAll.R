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

TYPE2=c('1142243F','1160920F','CID4290','CID4465','CID44971','CID4535')
DIR='./data/spatial_processed/'


OUT=.testAll_scRNA(TYPE1, TYPE2, DIR, .test_ikarus)
saveRDS(OUT, paste0(DIR,'/','eval_ikarus.rds'))




DIR='./data/spatial_processed/'
O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwo.rds'))
O3=readRDS(paste0(DIR,'/','eval_music.rds'))
O4=readRDS(paste0(DIR,'/','eval_ikarus.rds'))

OOO=cbind(O1[,4],O2[,4],O3[,4],O4[,4])
rownames(OOO)=paste0(O1[,1],'_',O1[,2])

pdf('plot/p13_boxplot_spatialRNA.pdf',width=3,height=3.5)
COL=c(rep('indianred2',1),'royalblue1','grey70','grey70')
boxplot(OOO, pch='+',col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
barplot(OOO['BRCA_CID4465',],col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
barplot(OOO['CESC_CID4465',],col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
dev.off()


apply(OOO,2,mean)
#[1] 0.7245690 0.7158861 0.6040725 0.6087251

print(OOO['BRCA_CID4465',])
#0.9279730 0.9250000 0.9242374 0.8697060


print(OOO['CESC_CID4465',])
#0.8911914 0.8812555 0.8159483 0.8050619

#plot(O1[,4],O2[,4])
#abline(a=0,b=1)

#t.test(O1[,4],O2[,4],paired=T,alternative='greater')







