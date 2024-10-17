import torch
import matplotlib.pyplot as plt
import numpy as np
from esm import ProteinBertModel, ProteinBertTokenizer

# Step 1: Install and import the required libraries
# You need to install esm: pip install fair-esm
# If not installed, you can use: pip install torch matplotlib

# Load ESM-3 model and tokenizer
model_name = "esm2_t33_650M_UR50S"  # Adjust this if using a different model
model, alphabet = ProteinBertModel.from_pretrained(model_name), ProteinBertTokenizer.from_pretrained(model_name)

# Function to generate embeddings for a protein sequence
def get_embeddings(sequence):
    # Encode sequence
    batch = alphabet.get_batch_converter()([(sequence,)])
    tokens = torch.tensor(batch[0][0][1:-1]).unsqueeze(0)  # Skip start and end tokens

    # Get embeddings
    with torch.no_grad():
        results = model(tokens)
    embeddings = results['representations'][0].squeeze(0).numpy()

    return embeddings

# Compute contact map from embeddings
def compute_contact_map(embeddings):
    # Compute pairwise distances
    distance_matrix = np.linalg.norm(embeddings[:, None] - embeddings[None, :], axis=2)
    # Convert distances to similarity
    contact_map = np.exp(-distance_matrix)
    return contact_map

# Plot contact map
def plot_contact_map(contact_map, sequence):
    plt.figure(figsize=(10, 8))
    plt.imshow(contact_map, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Similarity')
    plt.title(f'Contact Map for Sequence of Length {len(sequence)}')
    plt.xlabel('Residue Position')
    plt.ylabel('Residue Position')
    plt.show()

# Example usage
sequence = "MADEEKLPPG"  # Replace with your protein sequence
embeddings = get_embeddings(sequence)
contact_map = compute_contact_map(embeddings)
plot_contact_map(contact_map, sequence)
