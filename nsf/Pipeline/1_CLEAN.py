import pandas as pd 

raw_file = '/Users/donovanvincent/Desktop/nsf/Data/therasabdab.csv'
sasa_file = '/Users/donovanvincent/Desktop/nsf/Data/0_Masked_Therapeutic_Data.csv' 

raw_df = pd.read_csv(raw_file)
sasa_df = pd.read_csv(sasa_file)

# Only take the Therasabdab seqs that have identified names in sasa file
clean_raw_df = raw_df[raw_df['Therapeutic'].isin(sasa_df['Therapeutic'])]
print(len(clean_raw_df)) #Double check

# Select out the names and sequences only
df_of_interest=clean_raw_df[['Therapeutic', 'Heavy Sequence']]

# Add on the SASA scores from the other file
merged_df = pd.merge(df_of_interest, sasa_df[['Ground_Truth_SASAs']], left_index=True, right_index=True)

# Send to a new file 
merged_df.to_csv('/Users/donovanvincent/Desktop/nsf/Data/4_Therapeutic_SASA.csv', index=False)