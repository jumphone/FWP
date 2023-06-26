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

OUT=.testAll_scRNA(TYPE1, TYPE2, DIR, .test_fwp)
saveRDS(OUT, paste0(DIR,'/','eval_fwp.rds'))
OUT=.testAll_scRNA(TYPE1, TYPE2, DIR, .test_fwo)
saveRDS(OUT, paste0(DIR,'/','eval_fwo.rds'))
OUT=.testAll_scRNA(TYPE1, TYPE2, DIR, .test_music)
saveRDS(OUT, paste0(DIR,'/','eval_music.rds'))


DIR='./data/spatial_processed/'
O1=readRDS(paste0(DIR,'/','eval_fwp.rds'))
O2=readRDS(paste0(DIR,'/','eval_fwo.rds'))

wilcox.test(O1[,4],O2[,4],paired=T,alternative='greater')


O3=readRDS(paste0(DIR,'/','eval_music.rds'))





