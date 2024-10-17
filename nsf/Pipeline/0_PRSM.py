'''
NOTE: This file was ran on rockfish with previously generated and
relaxed TheraSabDab structures using IgFold. As a result, file and 
directory paths will be different if trying to reimplement
'''

import os
import pandas as pd
import pyrosetta
from pyrosetta import *

# Rosetta Code to calculate per residue sasa metric (prsm)
pyrosetta.init("-mute all", silent=True)
prsm = pyrosetta.rosetta.core.simple_metrics.per_residue_metrics.PerResidueSasaMetric()

# Read in file and define the directory with pdbs
file = '/home/donnie/scr16-jgray21/donnie/therasabdab.csv'
df = pd.read_csv(file)
directory = '/home/donnie/scr16-jgray21/donnie/moderna/TheraSabDab/TheraSabDab_IgFold_Predictions/pdbs'

names = []
all_prsm = [] #Will be a list of lists

for file in sorted(os.listdir(directory)):
    temp = os.path.join(directory,file) # Create the file name 
    name = file[:-4] # Create the therapeutic name
    temp_pose = pyrosetta.pose_from_pdb(temp) # make a pose
    
    # per residue sasa
    per_res_sasa = prsm.calculate(temp_pose)
    sasa_list = list(per_res_sasa.values()) # Take the values out of dict and add to list


    # Append to this list
    names.append(name)
    all_prsm.append(sasa_list)
    print(f"Finished with {name}")


    # Adding each one in order method 
    data = {
        "Therapeutic": name,
        "Ground_Truth_SASAs": [sasa_list] #Wrap in a final list to call 0 and get the list of list back

    }

    new_df = pd.DataFrame(data)

    file_name = 'GT_SASAs2.csv'
    header = os.path.exists(file_name)
    new_df.to_csv(file_name, index=False, header= not header, mode='a')


# Try the list at the end method
data2 = {
    "Therapeutic": names,
    "Ground_Truth_SASAs": all_prsm
}

temp_df = pd.DataFrame(data2)
temp_df.to_csv('GT_SASAs.csv', index=False)