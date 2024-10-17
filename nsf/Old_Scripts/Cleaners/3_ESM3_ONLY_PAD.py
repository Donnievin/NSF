#NOTE: Should all have a germline since they come from a human repertoire
import pandas as pd 
from abnumber import Chain #conda install -c bioconda abnumber

# esm3 pads are _ for each character
# esm2 pads are <mask> for entire region

file = '/Users/donovanvincent/Desktop/nsf/Generation/therasabdab.csv'
df = pd.read_csv(file)
clean_df = df[df['Heavy Sequence'] != 'na'] #Remove mess ups

new_df = clean_df[['Therapeutic', 'Heavy Sequence']]#Only want the seqs columns



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

cdr1_len, cdr2_len, cdr3_len, masked_seqs = [], [], [], []

for seq in new_df['Heavy Sequence']:

    # Identify the ONCE cdr3 using abnumber [could choose other regions]
    chain=Chain(seq, 'imgt')
    region1 = chain.cdr1_seq
    region2 = chain.cdr2_seq
    region3 = chain.cdr3_seq

    # Keep track of the len of the cdr to validate later w/ cdr infill
    cdr1_len.append(len(region1))
    cdr2_len.append(len(region2))
    cdr3_len.append(len(region3))

    # Mask using the 2 different methods 
    mask1 = esm3_pad(seq, region1)
    mask2 = esm3_pad(mask1, region2)
    mask3 = esm3_pad(mask2, region3)

    masked_seqs.append(mask3)
    
new_df['cdr1_len'] = cdr1_len
new_df['cdr2_len'] = cdr2_len
new_df['cdr3_len'] = cdr3_len
new_df['masked_seq'] = masked_seqs

new_df.to_csv('/home/donnie/scr16-jgray21/donnie/Projects/nsf/Data/Masked_Therapeutic_data.csv', index=False)