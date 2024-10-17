import pandas as pd 
import matplotlib.pyplot as plt 

'''
Notes: For each L x 512 tensor, EACH residue is represented by 512 floats (AntiBERTy final layer)
Command: python [file] [tSNE or UMAP] [Seq Col] [Label col]
Example: python tSNE_AntiBERTy.py 0 heavy_mini.csv Heavy_Sequence Temperature
'''

import torch
#import umap #pip install umap-learn
import numpy as np
from sklearn.manifold import TSNE # pip install -U scikit-learn
from antiberty import AntiBERTyRunner #pip install antiberty



#Embed the list sequences using AntiBERTy
def antiberty_embeddings(list_of_sequences):
    antiberty = AntiBERTyRunner()
    embeddings = antiberty.embed(list_of_sequences) #dtype = list of size nxlx512

    #MIGHT NOT NEED GENERAL DIMENSIONS ANYMORE SINCE NO RESHAPE
    n = len(list_of_sequences) # Num of sequences in the list
    l =  max(len(seq) for seq in list_of_sequences)+2+1 #Len of longest sequence + 2 (antiberty paper) + 1 (0-indexed)
    d = 512 #Output heads from antiberty

    return embeddings

#Reshape the embeddings by mean pooling
def mean_pool_embeddings(embeddings):
    mean_pooled_embeddings = [] #Empty list to hold the mean pooled embeddings
    for tens in embeddings: 
        mean_pooled_embeddings.append(tens.mean(dim=0)) #Calculate the avg and add to the list
    new_embeddings = torch.stack(mean_pooled_embeddings) #Turn the list into a tensor now that shapes match
    return new_embeddings

def tsne_graph(data, labels):
    tsne = TSNE(n_components=2) #Perplexity changes the std term, higher = more neighbors
    tsne_data = tsne.fit_transform(data)

    # Convert labels to a numeric format if they're not already
    #unique_labels = list(set(labels))
    unique_labels = labels.unique()
    label_to_color = {label: i for i, label in enumerate(unique_labels)}
    color_indices = np.array([label_to_color[label] for label in labels])

    # Define a colormap
    cmap = plt.get_cmap('viridis', len(unique_labels))

    # Plot t-SNE results
    plt.figure(figsize=(10, 8))
    plt.scatter(tsne_data[:, 0], tsne_data[:, 1], c=color_indices, cmap=cmap, marker='o')
    plt.title('Method 4: Mean Pooled Embeddings t-SNE')
    plt.xlabel('t-SNE Dimension 1')
    plt.ylabel('t-SNE Dimension 2')

    # Create a legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(i), markersize=10, label=label) for i, label in enumerate(unique_labels)]
    plt.legend(handles=handles, loc='best', shadow=False)
    plt.show()


    return 0

def UMAP_graph(data, labels):
    # Apply UMAP
    reducer = umap.UMAP(n_components=2)  # You can set n_components=3 for 3D visualization
    umap_data = reducer.fit_transform(data)

    # Convert labels to a numeric format if they're not already
    #unique_labels = list(set(labels))
    unique_labels = labels.unique()
    label_to_color = {label: i for i, label in enumerate(unique_labels)}
    color_indices = np.array([label_to_color[label] for label in labels])

    # Define a colormap
    cmap = plt.get_cmap('viridis', len(unique_labels))

    # Plot t-SNE results
    plt.figure(figsize=(10, 8))
    plt.scatter(umap_data[:, 0], umap_data[:, 1], c=color_indices, cmap=cmap, marker='o')
    plt.title('Mean Pooled Embeddings UMAP')
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')

    # Create a legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(i), markersize=10, label=label) for i, label in enumerate(unique_labels)]
    plt.legend(handles=handles, loc='best', shadow=False)
    plt.show()
    return 0 

def push(seqs, labels):
    #seqs = combined_df['Seq'].values.tolist()
    
    #Embed (use an if statment to include newer embedding methods)
    embeddings = antiberty_embeddings(seqs)

    #Mean pool
    mpe = mean_pool_embeddings(embeddings)

    #Graph
    tsne_graph(mpe, labels)
    #UMAP_graph(mpe, labels)


# Give files some names
esm2_file = '/Users/donovanvincent/Desktop/nsf/ESM_Scripts/esm2_generation.csv'
esm3_file = '/Users/donovanvincent/Desktop/nsf/Generation/esm3_generation_w_end.csv'
therasabdab_file = '/Users/donovanvincent/Desktop/nsf/Generation/therasabdab.csv'

# Read in with pandas, and add a column to therasabdab to match format
esm2_df = pd.read_csv(esm2_file)
esm3_df = pd.read_csv(esm3_file)
therasabdab_df = pd.read_csv(therasabdab_file)
therasabdab_df['Seq_tag'] = ['na' for _ in range(len(therasabdab_df))]

# Only take out the columns we care about
new_therasabdab_df = therasabdab_df[['Seq_tag','Heavy Sequence']]

# Rename for formatting issues 
new_therasabdab_df = new_therasabdab_df.rename(columns={'Heavy Sequence':'Seq'})

# Give each of these some labels
esm2_df['Model'] = 'ESM2'
esm3_df['Model'] = 'ESM3'
new_therasabdab_df['Model'] = 'Therasabdab'


#new_esm2_df = esm2_df.iloc[::2] #Starting at 1, then every other
new_esm3_df = esm3_df.drop(columns=['secondary_str','sasa','function','coordinates','plddt','ptm','SeqOCon'])

#print(new_esm2_df.columns, new_esm3_df.columns, new_therasabdab_df.columns)

combined_df = pd.concat([esm2_df, new_esm3_df, new_therasabdab_df])
#combined_df.to_csv('yes.csv', index=False)
push(combined_df['Seq'],combined_df['Model'])