
# Score calculation using FW-based PCA

library(data.table)
library(irlba)
library(pROC)


.simple_combine <- function(exp_sc_mat1, exp_sc_mat2, FILL=FALSE){
    FILL=FILL
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    ##############################################
    if(FILL==TRUE){
        gene1=rownames(exp_sc_mat)
        gene2=rownames(exp_ref_mat)
        gene12=gene2[which(!gene2 %in% gene1)]
        gene21=gene1[which(!gene1 %in% gene2)]
        exp_sc_mat_add=matrix(0,ncol=ncol(exp_sc_mat),nrow=length(gene12))
        rownames(exp_sc_mat_add)=gene12
        colnames(exp_sc_mat_add)=colnames(exp_sc_mat)
        exp_ref_mat_add=matrix(0,ncol=ncol(exp_ref_mat),nrow=length(gene21))
        rownames(exp_ref_mat_add)=gene21
        colnames(exp_ref_mat_add)=colnames(exp_ref_mat)
        exp_sc_mat=rbind(exp_sc_mat, exp_sc_mat_add)
        exp_ref_mat=rbind(exp_ref_mat, exp_ref_mat_add)
    }
    ############################################
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }


.loadFileNoGap <-function(input_path){
   library(data.table)
   HEADER=as.character(fread(input_path,header=FALSE, nrows=1))
   HEADER=HEADER[2:length(HEADER)]
   NET=fread(input_path,header=FALSE,sep='\t',select=c(1), skip=1)$V1
   SIGNAL=fread(input_path,header=FALSE,sep='\t',select=c(2:(length(HEADER)+1)), skip=1)
   SIGNAL=as.matrix(SIGNAL)
   rownames(SIGNAL)=NET
   colnames(SIGNAL)=HEADER
   return(SIGNAL)
   }

.loadFile <-function(input_path){
   library(data.table)
   HEADER=as.character(fread(input_path,header=FALSE, nrows=1))
   NET=fread(input_path,header=FALSE,sep='\t',select=c(1))$V1
   SIGNAL=fread(input_path,header=FALSE,sep='\t',select=c(2:(length(HEADER)+1)))
   SIGNAL=as.matrix(SIGNAL)
   rownames(SIGNAL)=NET
   colnames(SIGNAL)=HEADER
   return(SIGNAL)
   }


.denCenter <-function(x){
    x.den=density(x)
    y=x.den$x[order(-x.den$y)[1]]
    return(y)
    }


.calFW<-function(data, tag ){
    D1=as.matrix(data)
    TAG=tag
    options(warn=-1)
    COR=cor(t(D1), TAG)[,1]
    options(warn=1)
    COR[which(is.na(COR))]=0
    FW=COR
    FW=FW-.denCenter(FW)
    return(FW)
    }


.calBFW<-function(fw){
    BFW=fw
    BFW[which(fw>0)]=1
    BFW[which(fw<=0)]=0
    return(BFW)
    }


fwo<-function(data, fw){
    D2=data
    FW=fw
    #############################
    INTER=intersect(names(FW),rownames(D2))
    UD2=D2[INTER,]
    load.1=FW[INTER]
    ##################################################
    Y=cor(UD2,load.1)[,1]
    names(Y)=colnames(UD2)
    return(Y)
    }


fwp <- function(data, fw, n=10){
    TITLE='# Score calculation using FW-based PCA #'
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    print(TITLE)
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    print(Sys.time()) 
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    print('starting...')
    #############################
    D2=as.matrix(data)
    D2=D2[which(matrixStats::rowVars(D2)>0),]
    FW=fw
    n=n
    ############################
    INTER=intersect(names(FW),rownames(D2))
    UD2=D2[INTER,]
    load.1=FW[INTER]
    ###########################
    print('calculating original score...')
    oY=fwo(UD2, load.1)
    ###########################
    print('scaling data...')
    SUD2=UD2
    RMS=sqrt(rowSums(SUD2**2)/(ncol(SUD2)-1))
    SUD2=SUD2/RMS
    ################
    print('conducting FW-based PCA')
    ######################
    XUD2=SUD2 * load.1
    ################
    fit.pca.fwp=irlba::prcomp_irlba( t(XUD2), n=n, center=FALSE, scale. = FALSE)
    this_pca=fit.pca.fwp$x
    vec=this_pca
    ##########################
    print('calculating predicted score...')
    vec.w=cor(vec, oY)[,1]
    pY.index=order(-abs(vec.w))[1]
    pY.d=vec.w[pY.index]/abs(vec.w[pY.index])
    pY=vec[,pY.index] * pY.d
    names(pY)=colnames(UD2)
    ###################################
    bestAbsCor=abs(vec.w[pY.index])
    print( paste0('bestAbsCor: ', round(bestAbsCor,2) ) )
    print( paste0('Index of bestAbsCor: ', pY.index) )
    ###################################
    load.2=.calFW(UD2,pY)
    ##################################
    print('calculating final score...')
    fit.pca.Y=prcomp(cbind(oY,pY),center=TRUE,scale.=TRUE)
    fY.o=fit.pca.Y$x[,1]
    fY.c=cor(fY.o, oY)
    fY.d=fY.c/abs(fY.c)
    fY=fY.o * fY.d
    ###################################
    print('finished!')
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    print(Sys.time())
    print(paste0(rep('#',nchar(TITLE)),collapse=''))
    ##################################
    return(fY)
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




