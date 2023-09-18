
<img src="https://fzhang.bioinfo-lab.com/img/tools/logo_fwp.png" height="250">

# Feature Weight Pro (FWP)

**Feature-weight based measurement of cancerous transcriptome using cohort-wide and sample-specific information**, **Cellular Oncology**, 2023, in press

# Content:

* [Usage](#usage)
* [Benchmark Dataset](#benchmark-dataset)

# Usage:

    # R version 4.2.0
    
    library('data.table') # 1.14.6
    library('irlba') # 2.3.5.1
    library('pROC') # 1.18.2
    
    source('fwp.R')
    
    TRAIN.TYPE='PCPG'
    TEST.TYPE='PCPG'
    
    TRAIN.DATA=.loadFileNoGap(paste0('./data/first/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0('./data/first/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    TEST.DATA=.loadFileNoGap(paste0('./data/first/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0('./data/first/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    
    # DATA: a matrix (Row: genes; Column: samples); after normalization (e.g. fpkm)
    # TAG: a vector (Tumor: 1; Normal: 0)
    
    # step1. calculate feature-weight (FW)
    FW=.calFW(TRAIN.DATA, TRAIN.TAG[,2])
    
    # OR, you can directly load our pre-trained FW file. 
    # https://github.com/jumphone/FWP/tree/main/feature_weight/TCGA_bulkRNA
    FW=readRDS('FW_bulk_TCGA_PCPG.rds') 
    
    # step2. calculate score
    out.fwp=fwp(TEST.DATA, FW)
    
    # step3. evaluation
    result=.evaluate(out.fwp, TEST.TAG[,2])   
    print(result)

# Benchmark Dataset:

### Pre-trained Feature-Weight Files:

RDS files: https://github.com/jumphone/FWP/tree/main/feature_weight/TCGA_bulkRNA

### Processed Data:

Baidu Cloud Storage：https://pan.baidu.com/s/1WxImyznSwDtlox7bh0EyGg?pwd=ilj2 

In the processed data, "correctDP" stands for "correct data phenotype".

In the "Baidu Cloud Storage", we provide the processed data of bulk & single-cell RNA-seq data.

The processed spatial data is too large to upload. 

### Raw Data: 

bulk RNA-seq (UCSC Xena): https://xenabrowser.net/datapages/

single-cell & spatial RNA-seq, Wu et al., Nat Genet, 2021: GSE176078 (GEO database)

single-cell RNA-seq, oligodendroglioma: GSE70630  

spatial RNA-seq, pancreatic ductal adenocarcinoma (PDAC): GSE211895 (raw counts); https://github.com/anvaly/SpatialPortal (metadata)

<img src="https://fzhang.bioinfo-lab.com/img/white.png" height="50">

-------------------------------------------------------------------------------------------------------------------

<img src="https://fzhang.bioinfo-lab.com/img/panda_happy_logo.png" height='150'>

#### More tools & studies: https://fzhang.bioinfo-lab.com/
