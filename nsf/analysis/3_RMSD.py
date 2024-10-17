import os
import numpy as np
import pandas as pd
from Bio import PDB
import analysis.old_analysis.rmsd as rmsd

def read_pdb_coordinates(filename):
    """Read PDB file and return a list of atomic coordinates."""
    coordinates = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract coordinates (columns 30-38, 38-46, 46-54)
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coordinates.append([x, y, z])
    return np.array(coordinates)

def get_heavy_chain(file):
    # Create a parser 
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('PDB_structure', file)

    HC_coords = []

    # Iterate through all chains in the structure
    for model in structure:
        for chain in model:
            # Check if the chain is a heavy chain
            if chain.id == 'H':
                for residue in chain:
                    for atom in residue:
                        HC_coords.append(atom.coord)
                        
    return HC_coords

def calculate_rmsd(coords1, coords2):
    """Calculate the RMSD between two sets of coordinates."""
    if coords1.shape != coords2.shape:
        raise ValueError("The two coordinate arrays must have the same shape.")
    
    return np.sqrt(np.mean(np.sum((coords1 - coords2) ** 2, axis=1)))

# Load the data file to save to 
file = '/Users/donovanvincent/Desktop/nsf/Data/3_All_Data_w_BPP.csv'
df = pd.read_csv(file)

# Lists to hold these values
no_sasa = []
sasa = []

# Steps: Read all file names in a directory
directory = '/Users/donovanvincent/Desktop/nsf/Data'
gt_dir = os.path.join(directory, 'pdbs')
no_SASA_dir = os.path.join(directory, 'No_SASA')
SASA_dir = os.path.join(directory, 'SASA')

for pdb in os.listdir(gt_dir):
    # Create full file names
    gt_pdb = os.path.join(gt_dir,pdb)
    no_SASA_pdb = os.path.join(no_SASA_dir,pdb)
    SASA_pdb = os.path.join(SASA_dir, pdb)

    try:
        # Read all the pdb coordinates
        gt_cords = get_heavy_chain(gt_pdb) # Special function b/c it has heavy and light chain
        no_sasa_cords = read_pdb_coordinates(no_SASA_pdb)
        sasa_cords = read_pdb_coordinates(SASA_pdb)
        
        # RMSD with no_sasa
        # no_sasa_rmsd = np.sqrt(np.mean(np.sum(gt_cords[:len(no_sasa_cords)] - no_sasa_cords)**2))
        no_sasa_rmsd = np.sqrt(np.mean(np.sum((gt_cords[:len(no_sasa_cords)]- no_sasa_cords) ** 2, axis=1)))
        no_sasa.append(no_sasa_rmsd)

        # RMSD with sasa 
        sasa_rmsd = np.sqrt(np.mean(np.sum((gt_cords[:len(sasa_cords)]- sasa_cords) ** 2, axis=1)))
        sasa.append(sasa_rmsd)

        if pdb == 'Alirocumab.pdb': 
            print(f"Alirocumab with no sasa has RMSD: {no_sasa_rmsd}")
            print(f"Alirocumab with sasa has RMSD: {sasa_rmsd}")

        #if pdb == ''


    except Exception as e:
        no_sasa.append("na")
        sasa.append("na")


data = {
    "no_SASA_RMSD": no_sasa,
    "SASA_RMSD": sasa
}

#new_df = pd.DataFrame(data)
#new_df.to_csv('/Users/donovanvincent/Desktop/nsf/Data/4_All_Data_w_RMSD.csv', index=False)

#df.to_csv('/Users/donovanvincent/Desktop/nsf/Data/4_All_Data_w_RMSD.csv', index=False)