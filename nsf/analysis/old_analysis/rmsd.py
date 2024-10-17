import numpy as np

def read_pdb(filename):
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

def calculate_rmsd(coords1, coords2):
    """Calculate the RMSD between two sets of coordinates."""
    if coords1.shape != coords2.shape:
        raise ValueError("The two coordinate arrays must have the same shape.")
    
    return np.sqrt(np.mean(np.sum((coords1 - coords2) ** 2, axis=1)))

def main(pdb1, pdb2):
    coords1 = read_pdb(pdb1)
    coords2 = read_pdb(pdb2)

    print(len(coords1), len(coords2))

    if len(coords1) != len(coords2):
        print("Warning: The number of atoms in the two PDB files does not match.")
    
    # Trimming to the smaller size
    min_length = min(len(coords1), len(coords2))
    rmsd_value = calculate_rmsd(coords1[:min_length], coords2[:min_length])

    print(f"RMSD: {rmsd_value:.3f} Å")

if __name__ == "__main__":
    pdb_file1 = "/Users/donovanvincent/Desktop/nsf/Data/pdbs/Mecbotamab.pdb"  # Replace with your first PDB file
    pdb_file2 = "/Users/donovanvincent/Desktop/nsf/Data/SASA/Mecbotamab.pdb"  # Replace with your second PDB file
    main(pdb_file1, pdb_file2)
