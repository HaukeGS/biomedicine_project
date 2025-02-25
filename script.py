import anndata as ad
import numpy as np
import scanpy as sc
from get_embedding import main_gene_selection
import pandas as pd
# Load the h5ad file
adata = ad.read_h5ad("/bigwork/nhwast8a/dataset/COVID19_ALL.h5ad")
# Load the scFoundation pre-trained model


# Generate embeddings
gene_list_df = pd.read_csv('../OS_scRNA_gene_index.19264.tsv', header=0, delimiter='\t')

gene_list = list(gene_list_df['gene_name'])

print("dir of file.h5ad")
print(dir(ad.read_h5ad("/bigwork/nhwast8a/dataset/COVID19_ALL.h5ad")))
print("list of genes of the model")

print(gene_list_df.head())
print("dir of the file main_gene_selecition")

print(dir(main_gene_selection))
print("var of data")

print(adata.var)
print("obsm of data")

print(adata.obsm)

print("Data")

print(adata)

embeddings,to_fill_columns, var = main_gene_selection(adata, gene_list)
np.save("embeddings.npy", embeddings)

adata.obsm["X_scfoundation"] = embeddings
# Compute neighbors and run UMAP
sc.pp.neighbors(adata, use_rep="X_scfoundation")
sc.tl.umap(adata)

# Plot UMAP
sc.pl.umap(adata, color=["cell_type"], title="scFoundation Embeddings")
