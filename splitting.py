import os
from tqdm import tqdm
import anndata
import numpy as np

# Große Datei laden (Backed-Modus für Speicheroptimierung)
h5ad_path = r"data/COVID19_ALL.h5ad"
adata = anndata.read_h5ad(h5ad_path, backed="r")

# Split-Einstellungen
num_splits = 10  # Anzahl der gewünschten Teil-Dateien
split_size = len(adata) // num_splits  # Anzahl der Zellen pro Datei

split_folder = r"data/dataset_split"
os.makedirs(split_folder, exist_ok=True)

# tqdm für Fortschrittsanzeige
for i in tqdm(range(num_splits), desc="Splitting", unit="file"):
    start = i * split_size
    end = (i + 1) * split_size if i < num_splits - 1 else len(adata)

    adata_sub = adata[start:end].to_memory()  # In den RAM laden
    split_path = os.path.join(split_folder, f"split_{i}.h5ad")
    adata_sub.write_h5ad(split_path)

print("✅ Splitting abgeschlossen!")