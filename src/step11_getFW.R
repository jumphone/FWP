

source('fwp.R')
DIR='./data/scRNA/'

TYPES=c(
       'BRCA','KIRP','LIHC','PCPG',
       'BLCA','CESC','CHOL','COAD',
       'ESCA','HNSC','KICH','LUSC',
       'PAAD','THYM','UCEC')

i=1
while(i<=length(TYPES)){
    TRAIN.TYPE=TYPES[i]
    print(TRAIN.TYPE)
    ####################
    TRAIN.DATA=.loadFileNoGap(paste0(DIR,'/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0(DIR,'/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    ####################
    FW=.calFW(TRAIN.DATA,TRAIN.TAG[,2])
    saveRDS(FW,paste0('./feature_weight/FW_bulk_TCGA_',TRAIN.TYPE,'.rds'))
    i=i+1
    }



