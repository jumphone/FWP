
<img src="https://fzhang.bioinfo-lab.com/img/tools/logo_fwp.png" height="250">

# Feature Weight Pro (FWP)

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
    
    # DATA: a matrix (Row: genes; Column: samples)
    # TAG: a vector (Tumor: 1; Normal: 0)
    
    # step1. calculate feature-weight (FW)
    FW=.calFW(TRAIN.DATA, TRAIN.TAG[,2])
    
    # step2. calculate score
    out.fwp=fwp(TEST.DATA, FW)
    
    # step3. evaluation
    result=.evaluate(out.fwp, TEST.TAG[,2])   
    print(result)
    
# Benchmark Dataset:

### Raw Data: 

bulk RNA-seq (UCSC Xena): https://xenabrowser.net/datapages/

scRNA-seq (Wu et al., Nat Genet, 2021): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078

### Processed Data:

Baidu Cloud Storage：https://pan.baidu.com/s/1WxImyznSwDtlox7bh0EyGg?pwd=ilj2 

In the processed data, "correctDP" stands for "correct data phenotype".

### Feature-Weight Files:

RDS files: https://github.com/jumphone/FWP/tree/main/feature_weight/TCGA_bulkRNA



<img src="https://fzhang.bioinfo-lab.com/img/white.png" height="50">

-------------------------------------------------------------------------------------------------------------------

<img src="https://fzhang.bioinfo-lab.com/img/panda_happy_logo.png" height='150'>

#### Explore more tools & studies: https://fzhang.bioinfo-lab.com/
