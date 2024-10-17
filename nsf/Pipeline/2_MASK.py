import pandas as pd 
from abnumber import Chain #conda install -c bioconda abnumber

# esm3 pads are _ for each character

file = '/Users/donovanvincent/Desktop/nsf/Data/0_Therapeutic_SASA.csv'
df = pd.read_csv(file)

def esm3_pad(seq, region):
    """
    Takes in a string, finds the CDR3, then replaces it with that many _ tokens for esm3
    """

    # Create a mask from the len of region
    mask = ''
    for aa in region: mask += '_' # esm3 pads are _ for each character

    # Swap current region with the mask 
    new_seq = seq.replace(region, mask) #replace the region in the str with the mask

    return new_seq 

fwr1_len, fwr2_len, fwr3_len, fwr4_len, masked_seqs = [], [], [], [], []

for seq in df['Heavy Sequence']:

    # Identify the ONCE cdr3 using abnumber [could choose other regions]
    chain=Chain(seq, 'imgt')
    region1 = chain.fr1_seq
    region2 = chain.fr2_seq
    region3 = chain.fr3_seq
    region4 = chain.fr4_seq

    # Keep track of the len of the cdr to validate later w/ cdr infill
    fwr1_len.append(len(region1))
    fwr2_len.append(len(region2))
    fwr3_len.append(len(region3))
    fwr4_len.append(len(region4))

    # Mask using the 2 different methods 
    mask1 = esm3_pad(seq, region1)
    mask2 = esm3_pad(mask1, region2)
    mask3 = esm3_pad(mask2, region3)
    mask4 = esm3_pad(mask3, region4)

    masked_seqs.append(mask4)

# Add all of the new columns to the df
df['fwr1_len'] = fwr1_len
df['fwr2_len'] = fwr2_len
df['fwr3_len'] = fwr3_len
df['fwr4_len'] = fwr4_len
df['masked_seq'] = masked_seqs

# Print this out to a new file
df.to_csv('/Users/donovanvincent/Desktop/nsf/Data/1_Therapeutic_SASA.csv', index=False)