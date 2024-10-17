#NOTE: Keep charge, instability index, 
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt 

# Study seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis #pip install biopython
#import Levenshtein # study LDs

file = '/Users/donovanvincent/Desktop/nsf/Data/2_All_Data.csv'
df = pd.read_csv(file)

def get_control():
    # Load file
    file = '/Users/donovanvincent/Desktop/nsf/Data/1_Fully_Prepped_Data.csv'
    temp_df = pd.read_csv(file)
    seqs = temp_df['masked_seq'].values.tolist()

    # Empty list to add after replacing
    fill_As = []

    # Fill _ w/ A's
    for i in seqs:
        fill_As.append(i.replace('_','G'))
    return fill_As

# Pick out all the sequences: Extra [] to keep it as a pd df
GT_seqs = df[['Heavy Sequence']]
no_sasa_seqs = df[['no_sasa_prediction']]
sasa_seqs = df[['w_sasa_prediction']]
control = get_control() # Fake predictions
control_seqs = pd.DataFrame({"Seq":control}) # core.frame.DF just like the other ones

# Rename for formatting issues 
GT_seqs = GT_seqs.rename(columns={'Heavy Sequence':'Seq'})
no_sasa_seqs = no_sasa_seqs.rename(columns={'no_sasa_prediction':'Seq'})
sasa_seqs = sasa_seqs.rename(columns={'w_sasa_prediction':'Seq'})

# Add a column to each for labels 
GT_seqs['Label'] = ['TheraSabDab' for i in range(len(GT_seqs))]
no_sasa_seqs['Label'] = ['No SASA' for i in range(len(no_sasa_seqs))]
sasa_seqs['Label'] = ['W SASA' for i in range(len(sasa_seqs))]
control_seqs['Label'] = ['Control (all G)' for i in range(len(control_seqs))]


# Combine the seqs
#combined_df = pd.concat([GT_seqs, no_sasa_seqs, sasa_seqs, control_seqs]) #
combined_df = pd.concat([no_sasa_seqs, sasa_seqs, GT_seqs,]) #

# lists for properties
a_list,b_list,c_list,d_list = [],[],[],[]

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


# Add all data to df 
combined_df['instability_index'] = a_list
combined_df['hydrophobicity'] = b_list
combined_df['isoelectric_pt'] = c_list
combined_df['charge_at_7'] = d_list
combined_df.to_csv('temp.csv',index=False)

# Cleaning 
cleaned_combined_df = combined_df[combined_df['hydrophobicity']!= 100]
cleaned_combined_df.to_csv('/Users/donovanvincent/Desktop/nsf/Data/3_All_Data_w_BPP.csv', index=False)

# Boxplot of Properties vs Model 
# Create a figure and a set of subplots
fig, axs = plt.subplots(2, 2, figsize=(15, 10))  # 2x2 grid of subplots

font = 28

# First boxplot
sns.boxplot(data=cleaned_combined_df, y='instability_index', hue='Label', ax=axs[0, 0])
axs[0, 0].set_title('Instability Index', fontsize=font)

# Second boxplot
sns.boxplot(data=cleaned_combined_df, y='hydrophobicity', hue='Label', ax=axs[0, 1])
axs[0, 1].set_title('Hydrophobicity', fontsize=font)

# Third boxplot
sns.boxplot(data=cleaned_combined_df, y='isoelectric_pt', hue='Label', ax=axs[1, 0])
axs[1, 0].set_title('Isoelectric Point', fontsize=font)

# Fourth boxplot
sns.boxplot(data=cleaned_combined_df, y='charge_at_7', hue='Label', ax=axs[1, 1])
axs[1, 1].set_title('Charge at pH 7', fontsize=font)

# Adjust layout
plt.tight_layout()

# Show the plot
plt.savefig('Test1.png')
plt.show()
