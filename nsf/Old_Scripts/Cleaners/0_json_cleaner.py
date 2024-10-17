import pandas as pd
import sys
import os
import json

#Idk
def clean_csv(sequences, metadata):
    run, link, author, species, bsource, btype, longitudinal, disease, subject, age, vaccine, chain, isotype = [], [], [], [], [], [], [], [], [], [], [], [], []

    metadata_interest = {
        'Run': run,
        'Link': link,
        'Author': author,
        'Species': species,
        'BSource': bsource,
        'BType': btype,
        'Longitudinal': longitudinal,
        'Disease': disease,
        'Subject': subject,
        'Age': age,
        'Vaccine': vaccine,
        'Chain': chain,
        'Isotype': isotype
    }

    for key in metadata_interest:
        metadata_interest[key].append(len(sequences) * [metadata[key]])

    for key in metadata_interest:
        sequences[key] = metadata_interest[key][0]

    return sequences

raw_csv_gz = "/home/donnie/scr16-jgray21/donnie/Projects/nsf/SRR14611333_1_Heavy_IGHG.csv.gz"

base_name = os.path.basename(raw_csv_gz)
#print(base_name)
clean_csv_pwd = os.path.splitext(base_name)[0]  # Remove .csv
#print(clean_csv_pwd)
clean_csv_pwd = os.path.splitext(clean_csv_pwd)[0]  # Remove .gz
#print(clean_csv_pwd)


#Every file is a big csv with like 90 names, but the metadata col is 1st because it helps to wrap it behind OAS
metadata_col_names = pd.read_csv(raw_csv_gz, nrows=0).columns #This only reads the metadata descriptors in the 1st col, but it is of type Index obj
metadata = ','.join(metadata_col_names) #Takes the col names and turns it into a single string
meta_data = json.loads(metadata) # Turns the json formatted str into a dictionary

#Get the names of every column after the 1st one
all_col_names_df = pd.read_csv(raw_csv_gz, nrows=1) #The 2nd row is the real col names

# All column names
#sequence,locus,stop_codon,vj_in_frame,v_frameshift,productive,rev_comp,complete_vdj,v_call,d_call,j_call,sequence_alignment,germline_alignment,sequence_alignment_aa,germline_alignment_aa,v_alignment_start,v_alignment_end,d_alignment_start,d_alignment_end,j_alignment_start,j_alignment_end,v_sequence_alignment,v_sequence_alignment_aa,v_germline_alignment,v_germline_alignment_aa,d_sequence_alignment,d_sequence_alignment_aa,d_germline_alignment,d_germline_alignment_aa,j_sequence_alignment,j_sequence_alignment_aa,j_germline_alignment,j_germline_alignment_aa,fwr1,fwr1_aa,cdr1,cdr1_aa,fwr2,fwr2_aa,cdr2,cdr2_aa,fwr3,fwr3_aa,fwr4,fwr4_aa,cdr3,cdr3_aa,junction,junction_length,junction_aa,junction_aa_length,v_score,d_score,j_score,v_cigar,d_cigar,j_cigar,v_support,d_support,j_support,v_identity,d_identity,j_identity,v_sequence_start,v_sequence_end,v_germline_start,v_germline_end,d_sequence_start,d_sequence_end,d_germline_start,d_germline_end,j_sequence_start,j_sequence_end,j_germline_start,j_germline_end,fwr1_start,fwr1_end,cdr1_start,cdr1_end,fwr2_start,fwr2_end,cdr2_start,cdr2_end,fwr3_start,fwr3_end,fwr4_start,fwr4_end,cdr3_start,cdr3_end,np1,np1_length,np2,np2_length,c_region,Redundancy,ANARCI_numbering,ANARCI_status
sequences_df = pd.read_csv(raw_csv_gz, header=1, usecols=['sequence_alignment_aa'])

sequences_df.to_csv('Master_Data.csv', index=False)

