#  Score calculation using FW-based PCA

    
    library('data.table')
    library('irlba')
    library('pROC')
    
    source('fwp.R')
    
    TRAIN.TYPE='PCPG'
    TEST.TYPE='PCPG'
    
    TRAIN.DATA=.loadFileNoGap(paste0('./data/56/Fitdevo/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0('./data/56/Fitdevo/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    TEST.DATA=.loadFileNoGap(paste0('./data/56/Fitdevo/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0('./data/56/Fitdevo/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    
    # Calculate feature-weight (FW)
    FW=.calFW(TRAIN.DATA, TRAIN.TAG[,2])
    
    # Calculate score
    out.fwp=fwp(TEST.DATA, FW)
    
    # Evaluation
    result=.evaluate(out.fwp, TEST.TAG[,2])
    
    print(result)
    
