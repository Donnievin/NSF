#NOTE: Should all have a germline since they come from a human repertoire
import pandas as pd 
from abnumber import Chain #conda install -c bioconda abnumber

# esm3 pads are _ for each character
# esm2 pads are <mask> for entire region

file = '/home/donnie/scr16-jgray21/donnie/Projects/nsf/Data/Data.csv'
df = pd.read_csv(file)


def esm2_pad(seq, region):
    """
    Takes in a string, finds the CDR3, then replaces that whole area with a '<mask>'
    """

    # Need a temp CHARACTER to mask the entire region or else the mask token will be spaced as well 
    temp = "?"

    # Swap current region with temp char [not in AA vocabulary]
    masked_seq = seq.replace(region, temp) #replace the region in the str with the temp

    # For more than 1 residue being masked, need to put a space btwn all letters 
    spaced_seq = ' '.join(masked_seq)

    # Swap that temporary mask with the correct mask needed by esm2 now that all letters are spaced 
    new_seq = spaced_seq.replace(temp, '<mask>') #replace temp with actual mask

    return new_seq

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

cdr_len, esm2_seqs, esm3_seqs = [], [], []

for seq in df['sequences']:

    # Identify the ONCE cdr3 using abnumber [could choose other regions]
    chain=Chain(seq, 'imgt')
    region = chain.cdr3_seq

    # Keep track of the len of the cdr to validate later w/ cdr infill
    cdr_len.append(len(region))

    # Mask using the 2 different methods 
    esm2_seqs.append(esm2_pad(seq, region))
    esm3_seqs.append(esm3_pad(seq, region))

df['cdr_len'] = cdr_len
df['esm2_seq'] = esm2_seqs
df['esm3_seq'] = esm3_seqs

df.to_csv('/home/donnie/scr16-jgray21/donnie/Projects/nsf/Data/Masked_data.csv', index=False)