# /home/toolkit/local/bin/python3

import urllib.request
import anndata
import pandas as pd
from pathlib import Path
from ikarus import classifier, utils, data, gene_list
import upsetplot
import scanpy
from fwp_source import *
####################
this_path='./data/revise_scRNA_glioma/'
TYPE=['BRCA','KIRP','LIHC','PCPG',
       'BLCA','CESC','CHOL','COAD',
       'ESCA','HNSC','KICH','LUSC',
       'PAAD','THYM','UCEC',
       'OLIG']
#####################
i=0
while i<len(TYPE):
    this_type=TYPE[i]
    print(this_type)
    this_h5_path=this_path+'/'+this_type+'/train.h5ad'
    this_adata=loadH5(this_h5_path)
    this_dir=this_path+'/'+this_type+'/ikarus/'
    getSig(this_adata, this_dir)
    ikarus_train(this_dir+'/signatures.gmt', this_adata, this_dir)
    i=i+1


###########################################
i=0
while i<len(TYPE):
    j=0
    while j <len(TYPE):
        this_train=TYPE[i]
        this_test=TYPE[j]
        print(this_train)
        print(this_test)
        this_dir=this_path+'/'+this_train+'/ikarus/'
        this_sig_path=this_dir+'/signatures.gmt'
        this_model_path=this_dir+'/core_model.joblib'
        this_h5_path=this_path+'/'+this_test+'/test.h5ad'
        this_adata=loadH5(this_h5_path)
        ikarus_pred(this_sig_path, this_model_path, this_adata, this_dir, this_test)
        j=j+1
    i=i+1






