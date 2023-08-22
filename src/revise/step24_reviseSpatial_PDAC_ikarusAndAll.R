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

TYPE2=c('PDAC_DA_FFPE_PR','PDAC_DB_FFPE_PR','PDAC_DC_FFPE_PR')
DIR='./data/revise_spatial_PDAC/'


OUT=.testAll_scRNA(TYPE1, TYPE2, DIR, .test_ikarus)
saveRDS(OUT, paste0(DIR,'/','eval_ikarus.rds'))




O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwo.rds'))
O3=readRDS(paste0(DIR,'/','eval_music.rds'))
O4=readRDS(paste0(DIR,'/','eval_ikarus.rds'))

OOO=cbind(O1[,4],O2[,4],O3[,4],O4[,4])
rownames(OOO)=O1[,1]

pdf('revisePlot/p01_boxplot_spatial.pdf',width=3,height=3.5)
COL=c(rep('indianred2',1),'royalblue1','grey70','grey70')
boxplot(OOO, pch='+',col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
barplot(OOO['BRCA',],col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
barplot(OOO['CESC',],col=COL,ylim=c(0,1))
abline(h=0.5,lty=2)
dev.off()


print(apply(OOO,2,mean))
#0.8443634 0.8436352 0.8078344 0.5760747





