#NOTE: Keep charge, instability index, 
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import Levenshtein

pre_BPP_df = pd.read_csv('/Users/donovanvincent/Desktop/nsf/Data/13_All_Data_w_BPP.csv')
BPP_df = pre_BPP_df[pre_BPP_df['isoelectric_pt']!='na']

pre_RMSD_df = pd.read_csv('/Users/donovanvincent/Desktop/nsf/Data/4_All_Data_w_RMSD.csv')
RMSD_df = pre_RMSD_df[pre_RMSD_df['no_SASA_RMSD']!= 'na']

# Reformat the RMSD, then rename and concat to make easier plot
no_SASA_df = RMSD_df[['no_SASA_RMSD']]
SASA_df = RMSD_df[['SASA_RMSD']]

no_SASA_df = no_SASA_df.rename(columns={'no_SASA_RMSD':'RMSD'})
SASA_df = SASA_df.rename(columns={'SASA_RMSD':'RMSD'})

no_SASA_df['Label'] = 'no_SASA'
SASA_df['Label'] = 'SASA'

reordered_RMSD = pd.concat([no_SASA_df, SASA_df])

# Calculate the LD for the final graph 
LD_df = pd.read_csv('/Users/donovanvincent/Desktop/nsf/Data/2_All_Data.csv')
no_sasa_LDs = [Levenshtein.distance(a,b) for a,b in zip(LD_df['Heavy Sequence'], LD_df['no_sasa_prediction'])]
labels = ['No Sasa' for i in no_sasa_LDs]

all_LDs = no_sasa_LDs + [Levenshtein.distance(a,b) for a,b in zip(LD_df['Heavy Sequence'], LD_df['w_sasa_prediction'])]
all_labels = labels + ["SASA" for i in range(len(all_LDs) - len(no_sasa_LDs))]

data = {
    'Seq Sim': all_LDs,
    'Label': all_labels
}

seq_sim_df = pd.DataFrame(data)



# Boxplot of Properties vs Model 
# Create a figure and a set of subplots
fig, axs = plt.subplots(2, 2, figsize=(15, 10))  # 2x2 grid of subplots

font = 28


# First boxplot
iep_values = BPP_df['isoelectric_pt'].astype(float).values.tolist()
sns.boxplot(data=BPP_df, y=iep_values, hue='Label', ax=axs[0, 0])
axs[0, 0].set_title('Instability Index', fontsize=font)
axs[0,0].legend(fontsize=font, loc='upper right')
axs[0,0].set_ylim(0,20)
axs[0,0].tick_params(axis='y', labelsize=15)

# Second boxplot
q_values = BPP_df['charge_at_7'].astype(float).values.tolist()
sns.boxplot(data=BPP_df, y=q_values, hue='Label', ax=axs[0, 1])
axs[0, 1].set_title('Charge at pH 7', fontsize=font)
axs[0,1].legend(fontsize=font, loc='upper right')
axs[0,1].set_ylim(-20,40)
axs[0,1].tick_params(axis='y', labelsize=15)

# Third boxplot
rmsd_values = reordered_RMSD['RMSD'].astype(float).values.tolist()
sns.boxplot(data=reordered_RMSD, y=rmsd_values, hue='Label', ax=axs[1, 0])
axs[1, 0].set_title('RMSD from PDB [Å]', fontsize=font)
axs[1,0].legend(fontsize=font, loc='upper right')
axs[1,0].tick_params(axis='y', labelsize=15)


# Fourth boxplot
sns.boxplot(data=seq_sim_df, y='Seq Sim', hue='Label', ax=axs[1, 1])
axs[1, 1].set_title('Number of Sequence Mutations', fontsize=font)
axs[1,1].legend(fontsize=font, loc='upper right')
axs[1,1].set_ylim(0,120)
axs[1,1].tick_params(axis='y', labelsize=15)


# Adjust layout
plt.tight_layout()

# Show the plot
plt.savefig('Test1.png')
plt.show()