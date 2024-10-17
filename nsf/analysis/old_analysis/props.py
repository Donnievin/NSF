#NOTE: Keep charge, instability index, 
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt 

# Study seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis #pip install biopython
from abnumber import Chain # To identify GL
import Levenshtein # study LDs

def align_truncate(target, query): 
    aligner = Align.PairwiseAligner()
        
    # Modifications to penalize gap generation
    aligner.mode = 'local'
    aligner.open_gap_score = -10  # Higher penalty for opening gaps
    aligner.extend_gap_score = -0.5  # Lower penalty for extending gaps

    alignments = aligner.align(target, query)
    best_alignment = alignments[0]
    tas = best_alignment.aligned[0][0][0] # Target alignment start
    tae = best_alignment.aligned[0][-1][-1] # Target alignment end
    qas = best_alignment.aligned[1][0][0] # query alignment start
    qae = best_alignment.aligned[1][-1][-1] # query alignment end
    aligned_target = target[tas:tae]
    aligned_query = query[qas:qae]

    return aligned_target, aligned_query


file = '/Users/donovanvincent/Desktop/nsf/Data/6_All_Data.csv'
df = pd.read_csv(file)

# Pick out all the sequences: Extra [] to keep it as a pd df
GT_seqs = df[['Heavy Sequence']]
no_sasa_seqs = df[['no_sasa_prediction']]
sasa_seqs = df[['w_sasa_prediction']]

# Rename for formatting issues 
GT_seqs = GT_seqs.rename(columns={'Heavy Sequence':'Seq'})
no_sasa_seqs = no_sasa_seqs.rename(columns={'no_sasa_prediction':'Seq'})
sasa_seqs = sasa_seqs.rename(columns={'w_sasa_prediction':'Seq'})

# Add a column to each for labels 
GT_seqs['Label'] = ['Ground Truth' for i in range(len(GT_seqs))]
no_sasa_seqs['Label'] = ['No SASA' for i in range(len(no_sasa_seqs))]
sasa_seqs['Label'] = ['W SASA' for i in range(len(sasa_seqs))]

# Combine the seqs
combined_df = pd.concat([no_sasa_seqs, GT_seqs, sasa_seqs])

# lists for properties
a_list,b_list,c_list,d_list, e_list, f_list = [],[],[],[],[],[]

for i in combined_df['Seq']:
    try:
        X = ProteinAnalysis(i)
        ii = X.instability_index()
        hphob = X.gravy()
        iep = X.isoelectric_point() #Maybe swap to molar_extinction_coefficient or secondary_structure_fraction
        neut_ch = X.charge_at_pH(7)

        a_list.append(ii)
        b_list.append(hphob)
        c_list.append(iep)
        d_list.append(neut_ch)

    except Exception as e:
        a_list.append(100)
        b_list.append(100)
        c_list.append(100)
        d_list.append(100)

    try: 
        chain = Chain(i, 'imgt')

        # Identify v and j genes using AbNumber
        vgene = str(chain.find_human_germlines(1)[0][0])
        jgene = str(chain.find_human_germlines(1)[1][0])

        # Calculate LDv
        v_target, v_query = align_truncate(vgene, i)
        LDv = Levenshtein.distance(v_target, v_query)
        e_list.append(LDv)

        # Calculate LDj
        j_target, j_query = align_truncate(jgene, i)
        LDj = Levenshtein.distance(j_target, j_query)
        f_list.append(LDj)
    
    except Exception as e: 
        e_list.append(100)
        f_list.append(100)

# Add all data to df 
combined_df['instability_index'] = a_list
combined_df['hydrophobicity'] = b_list
combined_df['isoelectric_pt'] = c_list
combined_df['charge_at_7'] = d_list
combined_df['LDvs'] = e_list
combined_df['LDjs'] = f_list

# Cleaning 
cleaned_combined_df = combined_df[(combined_df['hydrophobicity']!= 100) |(combined_df['LDjs']== 100)]

# Boxplot of Properties vs Model 
# Create a figure and a set of subplots
fig, axs = plt.subplots(2, 2, figsize=(15, 10))  # 2x2 grid of subplots

font = 28

# First boxplot
sns.boxplot(data=cleaned_combined_df, y='instability_index', hue='Label', ax=axs[0, 0])
axs[0, 0].set_title('Instability Index', fontsize=font)

# Second boxplot
sns.boxplot(data=cleaned_combined_df, y='charge_at_7', hue='Label', ax=axs[0, 1])
axs[0, 1].set_title('Charge at pH 7', fontsize=font)

# Third boxplot
sns.boxplot(data=cleaned_combined_df, y='LDvs', hue='Label', ax=axs[1, 0])
axs[1, 0].set_title('V Gene Edit Distance', fontsize=font)

# Fourth boxplot
sns.boxplot(data=cleaned_combined_df, y='LDjs', hue='Label', ax=axs[1, 1])
axs[1, 1].set_title('J Gene Edit Distance', fontsize=font)

# Adjust layout
plt.tight_layout()

# Show the plot
plt.savefig('Test1.png')
plt.show()