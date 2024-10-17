# token: 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 [DO NOT ERASE]
import os
import ast # Will be used to unpack the SASA values in the DF
import torch
import pandas as pd
from huggingface_hub import login
from esm.models.esm3 import ESM3 #pip install esm
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig

torch.manual_seed(31100) # For reproducibility

# Read in the prepped data file
prepped_file = '/Users/donovanvincent/Desktop/nsf/Data/1_Fully_Prepped_Data.csv'
df = pd.read_csv(prepped_file)

# Get the masked seqs and the SASA values for inference
masked_seqs = df['masked_seq'].values.tolist()
SASA_df = df['Ground_Truth_SASAs'].values #This is an array of strs (each str is 1 list)

# THANK GOD THIS WORKED, UGGGHHHHH
SASA_values = [ast.literal_eval(s) for s in SASA_df] # Neat package to turn it back to a list of list

# Will instruct you how to get an API key from huggingface hub, make one with "Read" permission.
#login()

# This will download the model weights and instantiate the model on your machine.
model: ESM3InferenceClient = ESM3.from_pretrained("esm3_sm_open_v1").to("cpu") # or "cpu"

no_sasa_recovered = []
sasa_recovered = [] 

for i, seq in enumerate(masked_seqs): 
    # Variable for the name to save pdbs 
    name = df['Therapeutic'][i]

    no_sasa_protein = ESMProtein(sequence=seq) #Only feed in sequence
    no_sasa_protein = model.generate(no_sasa_protein, GenerationConfig(track="sequence", num_steps=8, temperature=1))
    no_sasa_recovered.append(no_sasa_protein.sequence)
    
    # Predict the pdb
    no_sasa_protein = model.generate(no_sasa_protein, GenerationConfig(track="structure", num_steps=8))
    no_sasa_protein.to_pdb(f"/Users/donovanvincent/Desktop/nsf/Data/No_SASA/{name}.pdb")

    # Values do not need to be perfect, since they will be binned anyways
    sasa_protein = ESMProtein(sequence=seq, sasa=SASA_values[i][:len(seq)]) #Slap in the sasa values from IgFold pdbs predictions
    sasa_protein = model.generate(sasa_protein, GenerationConfig(track="sequence", num_steps=8, temperature=1))
    sasa_recovered.append(sasa_protein.sequence)
    # Predict the pdb
    sasa_protein = model.generate(sasa_protein, GenerationConfig(track="structure", num_steps=8))
    sasa_protein.to_pdb(f"/Users/donovanvincent/Desktop/nsf/Data/SASA/{name}.pdb")


df['no_sasa_prediction'] = no_sasa_recovered 
df['w_sasa_prediction'] = sasa_recovered 

df.to_csv('/Users/donovanvincent/Desktop/nsf/Data/2_All_Data.csv', index=False)
