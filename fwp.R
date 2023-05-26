
# Score calculation using FW & PCA

library(data.table)
library(HiClimR)
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


.generate_mean <- function(exp_sc_mat, TAG, print_step=500){
    exp_sc_mat=exp_sc_mat
    TAG=TAG

    print_step=print_step

    NewRef=matrix(0,ncol=length(unique(TAG)),nrow=nrow(exp_sc_mat))

    TAG=as.character(TAG)
    refnames=unique(TAG)
    total_num=length(refnames)
    outnames=c()
    i=1
    while(i<=length(refnames)){
        one=refnames[i]
        this_col=which(TAG==one)
        outnames=c(outnames,one)
        if(length(this_col) >1){
            this_new_ref=rowSums(exp_sc_mat[,this_col])/length(this_col)
            }else{
            this_new_ref = exp_sc_mat[,this_col]
            }
        NewRef[,i]=this_new_ref
        if(i%%print_step==1){print(paste0(i,' / ' ,total_num ))}
        i=i+1
        }
    print(paste0(total_num, ' / ', total_num ))
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }


.calTFIDF <- function(mat, scale.factor =10000) {
      mat <- as.matrix(x = mat)
      cSum <- colSums(x = mat)
      tf <- tcrossprod(x = mat, y = Matrix::Diagonal(x = 1 / cSum))
      rSum <- rowSums(x = mat)
      idf <- ncol(x = mat) / rSum
      nmat <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% tf
      nmat <- as.matrix( log1p( nmat  * scale.factor) ) 
      rownames(nmat)=rownames(mat)
      colnames(nmat)=colnames(mat)
      nmat[which(is.na(nmat))]=0
      return(nmat)
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



fwo<-function(data, fw){
    D2=data
    FW=fw
    #################################
    INTER=intersect(names(FW),rownames(D2))
    UD2=D2[INTER,]
    load.1=FW[INTER]
    ##################################################
    Y=cor(UD2,load.1)[,1]
    names(Y)=colnames(UD2)
    return(Y)
    }


.dwnClst <-function(data, n, npcs=10, topvar=2000, seed=123){
    # using clusters to downscale data 
    DATA=data
    n=n
    npcs=npcs
    topvar=min(topvar,nrow(DATA))
    seed=seed
    #############################
    print('dwnClst: variable features...')
    VAR=matrixStats::rowVars(DATA)
    USED=which(rank(-VAR)<=topvar)
    UDATA=DATA[USED,]
    print('dwnClst: pca...')
    fit.pca=irlba::prcomp_irlba( t(UDATA), n=npcs, center=TRUE, scale. = FALSE)
    PCA=fit.pca$x
    print('dwnClst: k-means...')
    set.seed(123)
    KM=kmeans(PCA, centers=n)
    CLST=KM$cluster
    DATA.CLST=.generate_mean(DATA, CLST )
    DATA.CLST=DATA.CLST[,order(as.numeric(colnames(DATA.CLST)))]
    #################
    OUT=list()
    OUT$clst=CLST
    OUT$data=DATA.CLST
    return(OUT)
    }



fwp <- function(data, fw, npcs=10, nmax=NULL){
    TITLE='# Score calculation using FW & PCA #'
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
    nmax=nmax
    ############################
    UD2=D2
    if( !is.null(nmax) ){
        if(ncol(D2)>nmax){
            print('sample-size is larger than nmax, downscale data...')
            DC.OUT=.dwnClst(D2,nmax)
            UD2=DC.OUT$data
            }
        }     
    #############################
    UD2.V=matrixStats::rowVars(UD2)
    UD2=UD2[which( UD2.V > 0 ),]
    #############################
    print('calculating original score...')
    oY=fwo(UD2, FW) 
    ############################
    print('calculating correlation...')
    ###########################
    UD2.NORM=.calTFIDF(UD2)
    UD2.cor=HiClimR::fastCor(UD2.NORM)
    #UD2.cor=cor(UD2.NORM)
    ############################
    print('conducting PCA...')
    load.1=oY
    load.1.scale=scale(load.1)[,1]
    ############################
    SUD2=UD2.cor+1
    SUD2.NORM=.calTFIDF(SUD2)
    ##########################
    XUD2=SUD2.NORM * load.1.scale
    fit.pca.fwp=irlba::prcomp_irlba( t(XUD2), n=npcs, center=FALSE, scale. = FALSE)
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
    print('calculating final score...')
    fit.pca.Y=prcomp(cbind(oY,pY),center=TRUE,scale.=TRUE)
    fY.o=fit.pca.Y$x[,1]
    fY.c=cor(fY.o, oY)
    fY.d=fY.c/abs(fY.c)
    fY=fY.o * fY.d
    ##################################
    if(!is.null(nmax)){
        if(ncol(D2)>nmax){
            CLST=DC.OUT$clst
            fY.all=rep(0,length(CLST))
            fY.all=fY[CLST]
            fY=fY.all
            }
        }
    ###################################
    names(fY)=colnames(D2)
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




