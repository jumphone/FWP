

library(reticulate)
use_python('/home/toolkit/local/bin/python3',required=TRUE)

########################
library(pROC)
library(MuSiC)
library(Biobase)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(anndata)

source('fwp.R')
source('fwp_source.R')

source('fwp_source.R')

this_path='./data/spatial_processed/'

TYPE=c('BRCA','KIRP','LIHC','PCPG',
       'BLCA','CESC','CHOL','COAD',
       'ESCA','HNSC','KICH','LUSC',
       'PAAD','THYM','UCEC',
       '1142243F','1160920F','CID4290','CID4465','CID44971','CID4535')
MODE=c('train','test')

USED_GENE=.getUsedGene(this_path,TYPE,MODE)

set.seed(123)
i=1
while(i<=length(TYPE)){
j=1
while(j<=length(MODE)){
#####################################
this_type=TYPE[i]
this_mode=MODE[j]
print(this_type)
print(this_mode)
.preAnn(this_path, this_type, this_mode, USED_GENE, CUT=10)
#####################################
j=j+1}
i=i+1}






