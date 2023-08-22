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

TYPE2=c('OLIG')
DIR='./data/revise_scRNA_glioma/'

OUT=.testAll_scRNA(TYPE1, TYPE2, DIR, .test_fwp)
saveRDS(OUT, paste0(DIR,'/','eval_fwp.rds'))
OUT=.testAll_scRNA(TYPE1, TYPE2, DIR, .test_fwo)
saveRDS(OUT, paste0(DIR,'/','eval_fwo.rds'))
OUT=.testAll_scRNA(TYPE1, TYPE2, DIR, .test_music)
saveRDS(OUT, paste0(DIR,'/','eval_music.rds'))


