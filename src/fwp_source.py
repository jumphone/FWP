# /home/toolkit/local/bin/python3

import urllib.request
import anndata
import pandas as pd
from pathlib import Path
from ikarus import classifier, utils, data, gene_list
import upsetplot
import scanpy


def loadH5(this_path):
    adata = anndata.read_h5ad(this_path)
    #adata=scanpy.read_h5ad(this_path)
    adata.var.rename(columns={'features':'gene_symbol'}, inplace=True)
    adata.var_names = adata.var['gene_symbol']
    adata.var_names_make_unique()
    adata.obs.head(2)
    adata.var.head(2)
    adata.to_df().head(2)
    adata.raw=adata
    print('finished')
    return(adata)

def getSig(adata, this_dir):
    adatas={}
    adatas['ikarus']=adata
    names=["ikarus"]
    obs_names = ["tumor"]
    label_upregs = ["Tumor"]
    label_downregs = ["Normal"]
    ##################
    signatures = gene_list.create_all(
        label_upregs_list=label_upregs,
        label_downregs_list=label_downregs,
        adatas_dict=adatas,
        names_list=names,
        obs_names_list=obs_names,
        integration_fun=utils.intersection_fun,
        top_x=300
        )
    #######################
    contents = upsetplot.from_contents(signatures)
    tumor_genes_intersection = contents.loc[(True)].values.ravel().tolist()
    tumor_genes = tumor_genes_intersection
    ######################
    tumor_genes_union = []
    for i in signatures.values():
        tumor_genes_union += i
    #######################
    tumor_genes_union = list(set(tumor_genes_union)) # unique genes
    tumor_genes_union = []
    for i in signatures.values():
        tumor_genes_union += i
    #######################
    tumor_genes_union = list(set(tumor_genes_union)) # unique genes
    #########################################
    adatas={}
    adatas['ikarus']=adata
    names=["ikarus"]
    obs_names = ["tumor"]
    label_upregs = ["Normal"]
    label_downregs = ["Tumor"]
    #########################################
    signatures = gene_list.create_all(
        label_upregs_list=label_upregs,
        label_downregs_list=label_downregs,
        adatas_dict=adatas,
        names_list=names,
        obs_names_list=obs_names,
        integration_fun=utils.intersection_fun,
        top_x=300
        )
    #########################################
    contents = upsetplot.from_contents(signatures)
    #########################################
    normal_genes_union = []
    for i in signatures.values():
        normal_genes_union += i
    #########################################
    normal_genes_union = list(set(normal_genes_union)) # unique genes
    overlap = list(set(tumor_genes_union) & set(normal_genes_union))
    normal_genes = list(set(normal_genes_union) - set(overlap))
    #########################################
    gene_list.save_gmt([normal_genes, tumor_genes], ["Normal", "Tumor"], out_dir=this_dir)





def ikarus_train(sig_path, adata, this_dir):
    signatures_path = Path(sig_path)
    model = classifier.Ikarus(signatures_gmt=signatures_path, out_dir=this_dir)
    train_adata_list = [adata]
    train_names_list = ["ikarus"]
    obs_columns_list = ["tumor"]
    model.fit(train_adata_list, train_names_list, obs_columns_list, save=True)
    
def ikarus_pred(sig_path, model_path, adata, this_dir, this_label):
    model = classifier.Ikarus(signatures_gmt=sig_path,out_dir=this_dir)
    model.load_core_model(model_path)
    predictions = model.predict(adata, this_label, save=True)


