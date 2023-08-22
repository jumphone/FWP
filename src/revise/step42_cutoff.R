



########################
library(pROC)
library(MuSiC)
library(Biobase)
library(Seurat)

source('fwp.R')
source('fwp_source.R')

TYPE=c('BRCA','KIRP','LIHC','PCPG')
DIR='./data/first/'


.testCut_fwp<-function(TRAIN.TYPE, TEST.TYPE,DIR){
    
    ####################
    TRAIN.DATA=.loadFileNoGap(paste0(DIR,'/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0(DIR,'/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    TEST.DATA=.loadFileNoGap(paste0(DIR,'/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0(DIR,'/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    ####################
    FW=.calFW(TRAIN.DATA, TRAIN.TAG[,2])
    ########################
    this_out=fwp(TEST.DATA, FW)
    ###################
    OUT1=this_out
    TAG=TEST.TAG[,2]
    COR=cor(OUT1, TAG, method='spearman')
     #################
     this_resp=TAG
     this_pred=OUT1
     this_roc=pROC::roc(response=this_resp, predictor=this_pred, quiet =TRUE,direction='<')
     this_auc=pROC::auc(this_roc)[1]
     this_best=coords(this_roc, "best", ret = "threshold")
     #################
     OUT=list()
     OUT$cor=COR
     OUT$auc=this_best
    return(OUT)
    }


OUT=.testAll(TYPE, DIR, .testCut_fwp)


saveRDS(OUT, file='BestCutoff.rds')





