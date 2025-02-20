import anndata as ad
import scfoundation as scf
import numpy as np
import scanpy as sc

# Load the h5ad file
adata = ad.read_h5ad("/bigwork/nhwast8a/dataset/COVID19_ALL.h5ad")
# Load the scFoundation pre-trained model
model = scf.load_model("human")  # Use "mouse" if the data is from mice
# Generate embeddings
embeddings = scf.get_embeddings(adata, model)
np.save("embeddings.npy", embeddings)

adata.obsm["X_scfoundation"] = embeddings
# Compute neighbors and run UMAP
sc.pp.neighbors(adata, use_rep="X_scfoundation")
sc.tl.umap(adata)

# Plot UMAP
sc.pl.umap(adata, color=["cell_type"], title="scFoundation Embeddings")