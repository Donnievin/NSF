import umap
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from sklearn.manifold import TSNE
from antiberty import AntiBERTyRunner #pip install antiberty

file = '/Users/donovanvincent/Desktop/nsf/Data/2_All_Data.csv'
df = pd.read_csv(file)

# Pick out all the sequences
GT_seqs = df['Heavy Sequence'].values.tolist()
no_sasa_seqs = df['no_sasa_prediction'].values.tolist()
sasa_seqs = df['w_sasa_prediction'].values.tolist()

#Antiberty embeddings
antiberty = AntiBERTyRunner()

gt_embeddings = antiberty.embed(GT_seqs)
ESM3_no_sasa_embeddings = antiberty.embed(no_sasa_seqs)
EMS3_w_sasa_embeddings = antiberty.embed(sasa_seqs)

# MPE 
gt_mpe = torch.stack([tensor.mean(dim=0) for tensor in gt_embeddings])
esm3_no_sasa_mpe = torch.stack([tensor.mean(dim=0) for tensor in ESM3_no_sasa_embeddings])
esm3_w_sasa_mpe = torch.stack([tensor.mean(dim=0) for tensor in EMS3_w_sasa_embeddings])

#tSNE (0) or UMAP (1)
tSNE_o_UMAP = 0
if tSNE_o_UMAP == 0: 
    dim_reduc = TSNE(n_components=2)
    name="tSNE"
elif tSNE_o_UMAP == 1: 
    dim_reduc = umap.UMAP(n_components=2)
    name="UMAP"

# Graph 
gt_red = dim_reduc.fit_transform(gt_mpe)
esm3_no_sasa_red = dim_reduc.fit_transform(esm3_no_sasa_mpe)
esm3_w_sasa_red = dim_reduc.fit_transform(esm3_w_sasa_mpe)

# Order to match the BP properties order
#plt.figure(figsize=(10,10))
plt.scatter(esm3_no_sasa_red[:, 0], esm3_no_sasa_red[:, 1], alpha=0.5, label='No SASA')
plt.scatter(gt_red[:, 0], gt_red[:, 1], alpha=0.5, label='TheraSabDab')
plt.scatter(esm3_w_sasa_red[:, 0], esm3_w_sasa_red[:, 1], alpha=0.5, label='SASA')
plt.legend()
plt.title('UMAP of Ab-PLM Embeddings')
plt.xlabel('Dimension 1')
plt.ylabel('Dimension 2')
plt.savefig(f'{name}.png')