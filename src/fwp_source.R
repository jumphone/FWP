
.runMusic<-function( TRAIN.DATA, TRAIN.TAG, TEST.DATA){
    library(MuSiC)
    library(Biobase)
    library(Seurat)
    DATA=TRAIN.DATA
    TAG=cbind(colnames(TRAIN.DATA), TRAIN.TAG)
    CUT=3
    this_normal_index=which(TAG[,2]==0)
    ###########################
    if(length(this_normal_index)<CUT){
        print('##############################')
        print('Lack normal sample, add extra!!!')
        print('##############################')
        eNum=CUT-length(this_normal_index)
        print(paste0('addNum = ', eNum))
        print('##############################')
        ###############################
        if(length(this_normal_index)>1){
            e_index=sample(this_normal_index,eNum,replace=TRUE)
            }else{
            e_index=rep(this_normal_index,eNum)
            }
        eDATA=DATA[,e_index]
        eTAG=TAG[e_index,]
        DATA=cbind(DATA,eDATA)
        TAG=rbind(TAG, eTAG)
        colnames(DATA)=paste0(colnames(DATA),'.',c(1:ncol(DATA)))
        TAG[,1]=colnames(DATA)
        }
    #########################
    TRAIN.DATA=DATA
    TRAIN.TAG=TAG[,2]
    ##########################
    MAT=as.matrix(TRAIN.DATA)
    ###############
    set.seed(123)
    PD= cbind(TRAIN.TAG, sample(c(1,2,3,4,5,6),ncol(TRAIN.DATA),replace=T))
    ################
    rownames(PD)=colnames(MAT)
    colnames(PD)=c('type','sample')
    PD=as.data.frame(PD)
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(list(counts=MAT),colData=DataFrame(type=PD$type, sample=PD$sample),,metadata=PD)
    ES.REF=sce
    ###################
    MAT=as.matrix(TEST.DATA)
    GS= ExpressionSet(assayData=MAT)
    PD= cbind(colnames(MAT),colnames(MAT))
    rownames(PD)=colnames(MAT)
    colnames(PD)=c('type1','type2')
    PD=as.data.frame(PD)
    PD <- new("AnnotatedDataFrame", data=PD)
    ES <- ExpressionSet(assayData=MAT,phenoData=PD)
    ES.BULK=ES
    ####################
    set.seed(123)
    Est.bulk = music_prop(bulk.mtx = exprs(ES.BULK), sc.sce = ES.REF, clusters = 'type',samples='sample')
    return(Est.bulk)
    }
   
.evaluate<-function(score, tag){
     OUT1=score
     TAG=tag
     COR=cor(OUT1, TAG, method='spearman')
     #################
     this_resp=TAG
     this_pred=OUT1
     this_roc=pROC::roc(response=this_resp, predictor=this_pred, quiet =TRUE,direction='<')
     this_auc=pROC::auc(this_roc)[1]
     OUT=list()
     OUT$cor=COR
     OUT$auc=this_auc
     return(OUT)
     }







.test_fwo<-function(TRAIN.TYPE, TEST.TYPE,DIR){
    ####################
    TRAIN.DATA=.loadFileNoGap(paste0(DIR,'/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0(DIR,'/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    TEST.DATA=.loadFileNoGap(paste0(DIR,'/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0(DIR,'/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    ####################
    FW=.calFW(TRAIN.DATA, TRAIN.TAG[,2])
    this_out=fwo(TEST.DATA, FW)
    this_eva=.evaluate(this_out, TEST.TAG[,2])
    return(this_eva)
    }






.test_fwp<-function(TRAIN.TYPE, TEST.TYPE,DIR,...){
    ####################
    TRAIN.DATA=.loadFileNoGap(paste0(DIR,'/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0(DIR,'/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    TEST.DATA=.loadFileNoGap(paste0(DIR,'/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0(DIR,'/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    ####################
    FW=.calFW(TRAIN.DATA, TRAIN.TAG[,2])
    ########################
    this_out=fwp(TEST.DATA, FW,...)
    ###################
    this_eva=.evaluate(this_out, TEST.TAG[,2])
    return(this_eva)
    }




.test_fwo_fw<-function(TRAIN.TYPE, TEST.TYPE,DIR, FW){
    ####################
    TEST.DATA=.loadFileNoGap(paste0(DIR,'/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0(DIR,'/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    ####################
    this_out=fwo(TEST.DATA, FW)
    this_eva=.evaluate(this_out, TEST.TAG[,2])
    return(this_eva)
    }


.test_fwp_fw<-function(TRAIN.TYPE, TEST.TYPE,DIR,FW){
    ####################
    TEST.DATA=.loadFileNoGap(paste0(DIR,'/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0(DIR,'/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    ####################
    this_out=fwp(TEST.DATA, FW)
    this_eva=.evaluate(this_out, TEST.TAG[,2])
    return(this_eva)
    }



.test_music<-function(TRAIN.TYPE, TEST.TYPE,DIR){
    ####################
    TRAIN.DATA=.loadFileNoGap(paste0(DIR,'/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0(DIR,'/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    TEST.DATA=.loadFileNoGap(paste0(DIR,'/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0(DIR,'/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    ####################
    TRAIN.RN=rownames(TRAIN.DATA)
    TEST.RN=rownames(TEST.DATA)
    TRAIN.DATA=TRAIN.DATA[unique(TRAIN.RN),]
    TEST.DATA=TEST.DATA[unique(TEST.RN),]
    ######################
    tmp.music=.runMusic(TRAIN.DATA, TRAIN.TAG[,2], TEST.DATA)
    this_out=tmp.music$Est.prop.weighted[,'1']
    this_eva=.evaluate(this_out, TEST.TAG[,2])
    return(this_eva)
    }

.test_ikarus<-function(TRAIN.TYPE, TEST.TYPE, DIR){
    ####################
    TEST.TAG=read.csv(paste0(DIR,'/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    this_csv_path=paste0(DIR, '/',TRAIN.TYPE,'/ikarus/',TEST.TYPE,'/prediction.csv')
    this_csv=read.csv(this_csv_path)
    this_out=this_csv[,'final_pred_proba_Tumor']
    this_eva=.evaluate(this_out, TEST.TAG[,2])
    return(this_eva)
    }







.testAll<-function( TYPE, DIR,  .test, ...){
    OUT=c()
    i=1
    while(i<=length(TYPE)){
    j=1
    while(j<=length(TYPE)){
    ##########################@@@@@@@
    TRAIN.TYPE=TYPE[i]
    TEST.TYPE=TYPE[j]
    print(TRAIN.TYPE)
    print(TEST.TYPE)
    this_out=.test(TRAIN.TYPE,TEST.TYPE, DIR, ...)
    #########
    this_tab=unlist(this_out)
    ##########
    this_x=data.frame(
       train=TRAIN.TYPE,
       test=TEST.TYPE,
       cor=this_tab[1],
       auc=this_tab[2]
       )
    OUT=rbind(OUT, this_x)
    ################################
    j=j+1}
    i=i+1}
    rownames(OUT)=paste0('id',c(1:nrow(OUT)))
    return(OUT)
    }

.testAll_scRNA<-function( TYPE1, TYPE2, DIR,  .test, ...){
    
    OUT=c()
    i=1
    while(i<=length(TYPE1)){
    j=1
    while(j<=length(TYPE2)){
    ##########################@@@@@@@
    TRAIN.TYPE=TYPE1[i]
    TEST.TYPE=TYPE2[j]
    print(TRAIN.TYPE)
    print(TEST.TYPE)
    this_out=.test(TRAIN.TYPE,TEST.TYPE, DIR, ...)
    #########
    this_tab=unlist(this_out)
    ##########
    this_x=data.frame(
       train=TRAIN.TYPE,
       test=TEST.TYPE,
       cor=this_tab[1],
       auc=this_tab[2]
       )
    OUT=rbind(OUT, this_x)
    ################################
    j=j+1}
    i=i+1}
    rownames(OUT)=paste0('id',c(1:nrow(OUT)))
    return(OUT)
    }


.stat<-function(TRAIN.TYPE, TEST.TYPE,DIR){
    ####################
    TRAIN.DATA=.loadFileNoGap(paste0(DIR,'/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0(DIR,'/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    TEST.DATA=.loadFileNoGap(paste0(DIR,'/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0(DIR,'/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    ####################
    n.train=ncol(TRAIN.DATA)
    normal.train=length(which(TRAIN.TAG[,2]==0))
    n.test=ncol(TEST.DATA)
    normal.test=length(which(TEST.TAG[,2]==0))
    this_eva=c(n.train,normal.train,n.test,normal.test)
    return(this_eva)
    }


.statAll<-function( TYPE, DIR ){
    OUT=c()
    i=1
    while(i<=length(TYPE)){
    j=1
    while(j<=length(TYPE)){
    ##########################@@@@@@@
    TRAIN.TYPE=TYPE[i]
    TEST.TYPE=TYPE[j]
    print(TRAIN.TYPE)
    print(TEST.TYPE)
    this_out=.stat(TRAIN.TYPE,TEST.TYPE, DIR)
    #########
    ##########
    this_x=data.frame(
       train=TRAIN.TYPE,
       test=TEST.TYPE,
       n.train=this_out[1],
       normal.train=this_out[2],
       n.test=this_out[3],
       normal.test=this_out[4]
       )
    OUT=rbind(OUT, this_x)
    ################################
    j=j+1}
    i=i+1}
    rownames(OUT)=paste0('id',c(1:nrow(OUT)))
    return(OUT)
    }















########################################



.getUsedGene<-function(this_path, TYPE, MODE){
####################
ALLGENE=c()
i=1
while(i<=length(TYPE)){
j=1
while(j<=length(MODE)){
this_type=TYPE[i]
this_mode=MODE[j]
print(this_type)
print(this_mode)
DATA=.loadFileNoGap(paste0(this_path,'/',this_type,'/mat_',this_mode,'.tsv'))
TAG=read.csv(paste0(this_path,'/',this_type,'/CorrectDP_',this_mode,'.csv'),header=F)
#DATA=DATA[which(rowSums(DATA)>0),]
ALLGENE=c(ALLGENE, rownames(DATA))
##################
j=j+1}
i=i+1}
USED_GENE=names(which(table(ALLGENE)==length(MODE)*length(TYPE)))
##############################
return(USED_GENE)
}




.preAnn<-function(this_path, this_type,this_mode,CUT=10, USED_GENE){

####################################
DATA=.loadFileNoGap(paste0(this_path,'/',this_type,'/mat_',this_mode,'.tsv'))
TAG=read.csv(paste0(this_path,'/',this_type,'/CorrectDP_',this_mode,'.csv'),header=F)
DATA=DATA[which(rownames(DATA) %in% USED_GENE),]
#DATA=DATA[which(rowSums(DATA)>0),]
#######################################
if(this_mode=='train'){
    this_normal_index=which(TAG[,2]==0)
    CUT=CUT
    if(length(this_normal_index)<CUT){
        print('##############################')
        print('Lack normal sample, add extra!!!')
        print('##############################')
        eNum=CUT-length(this_normal_index)
        print(paste0('addNum = ', eNum))
        print('##############################')
        ###############################
        if(length(this_normal_index)>1){
            e_index=sample(this_normal_index,eNum,replace=TRUE)
            }else{
            e_index=rep(this_normal_index,eNum)
            }
        eDATA=DATA[,e_index]
        eTAG=TAG[e_index,]
        DATA=cbind(DATA,eDATA)
        TAG=rbind(TAG, eTAG)
        colnames(DATA)=paste0(colnames(DATA),'.',c(1:ncol(DATA)))
        TAG[,1]=colnames(DATA)
        }
    }
#######################################
pbmc=CreateSeuratObject(counts=DATA,assay = "RNA")
pbmc$tumor01=TAG[,2]
pbmc$tumor=rep('Tumor',length(TAG[,2]))
pbmc$tumor[which(TAG[,2]==0)]='Normal'
SaveH5Seurat(pbmc, filename = paste0(this_path,'/',this_type,'/', this_mode,'.h5Seurat'),overwrite=TRUE)
Convert( paste0(this_path,'/',this_type,'/', this_mode,'.h5Seurat'), dest = "h5ad", overwrite=TRUE)

####################
}




fwp_dw <- function(data, fw, npcs=2, dw=1, mode='final'){
    TITLE='# Score calculation using Feature Weight Pro #'
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    print(TITLE)
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    print(Sys.time())
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    print('starting...')
    #############################
    D2=as.matrix(data)
    FW=fw
    npcs=npcs
    mode=mode
    dw=dw
    ###########################
    INTER=intersect(rownames(D2),names(FW))
    INTER=sample(INTER, length(INTER)*dw )
    D2=D2[INTER, ]
    FW=FW[INTER]
    ############################
    UD2=D2
    UD2.V=matrixStats::rowVars(UD2)
    UD2=UD2[which( UD2.V > 0 ),]
    FW=FW[rownames(UD2)]
    #############################
    print('calculating original score...')
    oY=fwo(UD2, FW)
    ############################################
    print('calculating score(pca) & score(afw)...')
    ###################################
    # Global Information: Score (PCA)
    FIT.PCA=irlba::prcomp_irlba(t( UD2 * FW ), n=npcs, center = TRUE, scale. = FALSE)
    PCA=FIT.PCA$x
    FIT.LM=lm(oY~PCA)
    pY=predict(FIT.LM)
    names(pY)=names(oY)
    if(mode %in% c('global','pca','g','p')){return(pY)}
    ####################################
    # Local Information: Score (AFW)
    ZMAT=apply(UD2, 2, .pcc_perturb, FW, only_pos=TRUE)
    TMP=apply(ZMAT, 1, max)
    adjFW=TMP * sign(FW)
    adjFW=adjFW[which(adjFW!=0)]
    aY=fwo(UD2, adjFW)
    if(mode %in% c('local','afw','l','a')){return(aY)}
    ####################################
    print('calculating final score...')
    X=cbind(pY, aY)
    FIT=lm(oY~X)
    fY=predict(FIT)
    names(fY)=names(oY)
    ###################################
    print('finished!')
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    print(Sys.time())
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    ##################################
    return(fY)
    }
 





.test_fwp_dw<-function(TRAIN.TYPE, TEST.TYPE,DIR,...){
    ####################
    TRAIN.DATA=.loadFileNoGap(paste0(DIR,'/',TRAIN.TYPE,'/mat_train.tsv'))
    TRAIN.TAG=read.csv(paste0(DIR,'/',TRAIN.TYPE,'/CorrectDP_train.csv'),header=F)
    TEST.DATA=.loadFileNoGap(paste0(DIR,'/',TEST.TYPE,'/mat_test.tsv'))
    TEST.TAG=read.csv(paste0(DIR,'/',TEST.TYPE,'/CorrectDP_test.csv'),header=F)
    ####################
    FW=.calFW(TRAIN.DATA, TRAIN.TAG[,2])
    ########################
    this_out=fwp_dw(TEST.DATA, FW,...)
    ###################
    this_eva=.evaluate(this_out, TEST.TAG[,2])
    return(this_eva)
    }



