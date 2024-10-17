import umap
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from sklearn.manifold import TSNE
from antiberty import AntiBERTyRunner #pip install antiberty

file = '/Users/donovanvincent/Desktop/nsf/Data/6_All_Data.csv'
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

gt_mpe = gt_mpe.flatten().numpy()
esm3_no_sasa_mpe = esm3_no_sasa_mpe.flatten().numpy()
esm3_w_sasa_mpe = esm3_w_sasa_mpe.flatten().numpy()

# Histplot of the MPE
plt.hist(gt_mpe, color='red', alpha=0.5, bins=20, label='TheraSabDab')
plt.hist(esm3_no_sasa_mpe, color="blue", alpha=0.5, bins=20, label='No SASA')
plt.hist(esm3_w_sasa_mpe,color="green", alpha=0.5, bins=20, label='SASA')
plt.legend()
plt.title('Sequences from Different Tracks')
plt.xlabel('Mean Pooled Embeddings')
plt.ylabel('Frequency')
plt.savefig('hist.png')