
<img src="https://fzhang.bioinfo-lab.com/img/tools/logo_fwp.png" height="250">

# Score calculation using FWP (Feature Weight Pro)


# Usage:

    library('data.table')
    library('irlba')
    library('pROC')
    
    source('fwp.R')
    
    TRAIN.TYPE='PCPG'
    TEST.TYPE='PCPG'
    
    TRAIN.DATA=.loadFileNoGap(paste0('./data/first/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0('./data/first/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    TEST.DATA=.loadFileNoGap(paste0('./data/first/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0('./data/first/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    
    DATA: a matrix (Row: genes; Column: samples)
    TAG: a vector (Tumor: 1; Normal: 0)
    
    # step1. calculate feature-weight (FW)
    FW=.calFW(TRAIN.DATA, TRAIN.TAG[,2])
    
    # step2. calculate score
    out.fwp=fwp(TEST.DATA, FW)
    
    # step3. evaluation
    result=.evaluate(out.fwp, TEST.TAG[,2])   
    print(result)
    
