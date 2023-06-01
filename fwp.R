
# Score calculation using FW & PCA

library(data.table)
library(HiClimR)
library(irlba)
library(pROC)


.rmout<-function(x){
    this_q3=quantile(x,0.75)
    this_q1=quantile(x,0.25)
    this_up=this_q3+1.5*(this_q3-this_q1)
    this_lw=this_q1-1.5*(this_q3-this_q1)
    y=x
    y[which(x>this_up)]=this_up
    y[which(x<this_lw)]=this_lw
    return(y)
    }

.calABCD<-function(X, Y, X_base, Y_base){
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


.pcc_perturb_base0<-function(X, Y, only_pos=TRUE){
    only_pos=only_pos
    #############################
    X=.rmout(X)
    Y=.rmout(Y)
    N=length(X)
    ############################
    X_base = 0
    Y_base = 0
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


.pcc_perturb<-function(X, Y, only_pos=FALSE ){
    only_pos=only_pos
    #############################
    N=length(X)
    ############################
    X_base = mean(X)
    Y_base = mean(Y)
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


fwp <- function(data, fw, npcs=2){
    TITLE='# Score calculation using Feature-Weight Pro #'
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
    ###########################
    INTER=intersect(rownames(D2),names(FW))
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
    print('calculating predicted score...')
    ###############
    SUD2=apply(UD2,2,scale)
    ZMAT=apply(SUD2, 2, .pcc_perturb_base0, FW, only_pos=TRUE)
    TMP=apply(ZMAT, 1, max)
    adjFW=TMP * sign(FW)
    adjFW=adjFW[which(adjFW!=0)]
    pY=fwo(UD2, adjFW)
    ###################################
    FIT=irlba::prcomp_irlba(t( UD2 * FW ), n=npcs, center = TRUE, scale. = FALSE)
    PCA=FIT$x   
    cY=predict(lm(oY~PCA))
    ####################################
    print('calculating final score...')
    X=cbind(pY, cY)
    fY=predict(lm(oY~X))
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




