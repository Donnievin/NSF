import pandas as pd 

file = '/home/donnie/scr16-jgray21/donnie/Projects/nsf/Data/Master_Data.csv'

df = pd.read_csv(file)
seqs = df['sequence_alignment_aa']

new = []

for seq in seqs:
      if len(seq) > 130:
          new.append(seq)


new_df = pd.DataFrame(new)
new_df.to_csv("/home/donnie/scr16-jgray21/donnie/Projects/nsf/Data/Data.csv", index=False)