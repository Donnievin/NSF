import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt 

# Study seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis #pip install biopython

# Study GL
#import Levenshtein
#from abnumber import Chain #conda install -c bioconda abnumber


# Give files some names
esm2_file = '/Users/donovanvincent/Desktop/nsf/ESM_Scripts/esm2_generation.csv'
esm3_file = '/Users/donovanvincent/Desktop/nsf/Generation/esm3_generation_noend.csv'
therasabdab_file = '/Users/donovanvincent/Desktop/nsf/Generation/therasabdab.csv'

# Read in with pandas, and add a column to therasabdab to match format
esm2_df = pd.read_csv(esm2_file)
esm3_df = pd.read_csv(esm3_file)
therasabdab_df = pd.read_csv(therasabdab_file)
therasabdab_df['Seq_tag'] = ['na' for _ in range(len(therasabdab_df))]

# Only take out the columns we care about
new_therasabdab_df = therasabdab_df[['Seq_tag','Heavy Sequence']]

# Rename for formatting issues 
new_therasabdab_df = new_therasabdab_df.rename(columns={'Heavy Sequence':'Seq'})

# Give each of these some labels
esm2_df['Model'] = 'ESM2'
esm3_df['Model'] = 'ESM3'
new_therasabdab_df['Model'] = 'Therasabdab'


new_esm2_df = esm2_df.iloc[::2] #Starting at 1, then every other
new_esm3_df = esm3_df.drop(columns=['secondary_str','sasa','function','coordinates','plddt','ptm','SeqOCon'])

#print(new_esm2_df.columns, new_esm3_df.columns, new_therasabdab_df.columns)

combined_df = pd.concat([new_esm2_df, new_esm3_df, new_therasabdab_df])


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

# Cleaning 
cleaned_combined_df = combined_df[combined_df['hydrophobicity']!= 100]
#cleaned_combined_df.to_csv('test.csv', index=False)

# Boxplot of Properties vs Model 
# Create a figure and a set of subplots
fig, axs = plt.subplots(2, 2, figsize=(15, 10))  # 2x2 grid of subplots

font = 22

# First boxplot
sns.boxplot(data=cleaned_combined_df, y='instability_index', hue='Model', ax=axs[0, 0])
axs[0, 0].set_title('Instability Index', fontsize=font)

# Second boxplot
sns.boxplot(data=cleaned_combined_df, y='hydrophobicity', hue='Model', ax=axs[0, 1])
axs[0, 1].set_title('Hydrophobicity', fontsize=font)

# Third boxplot
sns.boxplot(data=cleaned_combined_df, y='isoelectric_pt', hue='Model', ax=axs[1, 0])
axs[1, 0].set_title('Isoelectric Point', fontsize=font)

# Fourth boxplot
sns.boxplot(data=cleaned_combined_df, y='charge_at_7', hue='Model', ax=axs[1, 1])
axs[1, 1].set_title('Charge at pH 7', fontsize=font)

# Adjust layout
plt.tight_layout()

# Show the plot
plt.savefig('Test1.png')
plt.show()