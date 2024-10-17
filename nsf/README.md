------------------------
To reproduce the results found in this proposal, run the following pipeline:
```bash
python Run_NSF_Pipeline.py
```

This will run all files in the order outlined in the pipline below: 


------------------------


##  NOTE: Must have a PyRosetta License (free to academics)
```bash
https://www.pyrosetta.org/home/licensing-pyrosetta
```

------------------------

# Pipeline
## TheraSabDab as Ground Truth (DONE)
1. Take TheraSabDab Heavy Sequences, measure per residue SASA using PyRosetta 
Code: ```Pipeline/0_PRSM.py```

## Clean up (DONE)
1. Take only the sequences that I was able to calculate the sasa for (a.k.a. Heavy Seq != na) 
Code: ```Pipeline/1_CLEAN.py```

## Masking (DONE)
1. Find and mask all 4 FWRs to prep the seqs for inference
Code:  ```Pipeline/2_MASK.py```

## ESM3 Inference, BPP vs no BPP (DONE)
1. Feed the each masked sequence through ESM3 twice, once without the SASA scores (as calculated above) and once with the SASA scores
Code: ```Pipeline/3_ESM3.py``` 

## Model Analysis (DONE)
1. UMAP - Feed all seqs through AntiBERTy, embed them and make a UMAP (tSNE did not separate as well) 
Code: ```analysis/1_Dim_Red.py```

2. BPPs of seqs - Feed all seqs through Seq based Oracles for comparison to therapeutic
Code:  ```analysis/2_Properties.py``` to see how they match up with the known therapeutics

3. (OPTIONAL) Use IgFold pRMSD as an indication for "how far from Ig-like" the structure is

4. (OPTIONAL) Compare the pdbs to their "ground truth" and calculate an RMSD between the two pdbs (DONE)
Code: ```analysis/3_RMSD.py```

5. Saved the data to a new_csv to make for faster analysis, remade graphs faster using new code
Code: ```analysis/9_Fast_Graphs.py```


# Properties 
Instability Index
Hydrophobicity
Isoelectric Pt (KEEP) [Median is the same]
Charge at pH7 
LDv + LDj => Distinct from germline

`Top Left -> Instability Index`
Top right -> 
Bottom Left -> 
`Bottom Right ->LDv + LDj`