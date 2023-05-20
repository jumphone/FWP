# Score calculation using FW-based PCA


library(data.table)
library(irlba)
library(pROC)

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

.rmOut<-function(X){
    X=X
    Q3=quantile(X,0.75)
    Q1=quantile(X,0.25)
    RANGE=Q3-Q1
    UP=Q3+1.5*RANGE
    LW=Q1-1.5*RANGE
    OUT=X
    OUT[which(X>UP)]=UP
    OUT[which(X<LW)]=LW
    return(OUT)
    }

.calABCD<-function(X, Y, X_base, Y_base){
    X=X
    Y=Y
    X_base=X_base
    Y_base=Y_base
    X_delta=X-X_base
    Y_delta=Y-Y_base
    A = sum(X_delta * Y_delta)
    B = sum(X_delta ** 2)
    C = sum(Y_delta ** 2)
    ############################
    D = A / sqrt( B * C )
    OUT=list()
    OUT[['A']]=A
    OUT[['B']]=B
    OUT[['C']]=C
    OUT[['D']]=D
    return(OUT)
    }

.pcc_perturb<-function(X, Y, only_pos=TRUE){
    X=.rmOut(X)
    Y=.rmOut(Y)
    only_pos=only_pos
    #############################
    N=length(X)
    ############################
    X_mean = mean(X)
    Y_mean = mean(Y)
    ###########################
    X_base = X_mean
    Y_base = Y_mean
    X_delta = X - X_base
    Y_delta = Y - Y_base
    ###########################
    ABCD = .calABCD(X, Y, X_base, Y_base)
    ###########################
    D = ABCD[['D']]
    D_plus = ( ABCD[['A']] + X_delta * Y_delta ) / sqrt( (ABCD[['B']] + X_delta**2) * (ABCD[['C']] + Y_delta**2) )
    ###########################
    M = D_plus - D
    S = (1-D**2)/(N-1)
    ###########################
    Z = M / S
    Z[which(is.na(Z))]=0
    if(only_pos==TRUE){Z[which(Z<0)]=0 }
    return( Z )
    }


.calFW<-function(data, tag){
    D1=data
    TAG=tag
    options(warn=-1)
    COR=cor(t(D1), TAG)[,1]
    options(warn=1)
    COR[which(is.na(COR))]=0
    FW=COR
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
    print('starting...')
    D2=as.matrix(data)
    FW=fw
    n=n
    #############################
    INTER=intersect(names(FW),rownames(D2))
    UD2=D2[INTER,]
    load.1=FW[INTER]
    ###########################
    print('calculating original Y...')
    oY=cor(UD2,load.1)[,1]
    ###########################
    print('scaling data...')
    SUD2=UD2
    RMS=sqrt(rowSums(UD2**2)/(ncol(UD2)-1))
    SUD2=SUD2/RMS
    SUD2[which(is.na(SUD2))]=0
    ################
    print('conducting FW-based PCA')
    XUD2=SUD2 * load.1
    ################
    fit.pca.2=irlba::prcomp_irlba( t(XUD2),
              n=n, center=FALSE, scale. = FALSE)
    vec=fit.pca.2$x
    ##########################
    print('calculating predicted Y...')
    load.2=cor(vec, oY)[,1]
    #print(load.2)
    pY.index=order(-abs(load.2))[1]
    pY.d=load.2[pY.index]/abs(load.2[pY.index])
    pY=vec[,pY.index] * pY.d
    names(pY)=colnames(UD2)
    ###################################
    bestAbsCor=abs(load.2[pY.index])
    print( paste0('bestAbsCor: ', round(bestAbsCor,2) ) )
    print( paste0('Index of bestAbsCor: ', pY.index) )
    ###################################
    options(warn=-1)
    load.3=cor(t(UD2),pY)[,1]
    options(warn=1)
    load.3[which(is.na(load.3))]=0
    ##################################
    print('calculating final Y...')
    Z=.pcc_perturb(load.1, load.3, only_pos=TRUE)
    load.4 = load.1 * Z
    fY=cor(UD2, load.4)[,1]
    ###################################
    print('finished!')
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













