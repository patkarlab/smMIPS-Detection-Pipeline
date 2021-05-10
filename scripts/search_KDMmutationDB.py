import pandas as pd
import os, sys

args = sys.argv

mutDb_file='/scratch/hematopath/hemoseq_v2/mutation_detector_nextflow/kdm_mutation_database/kdm_mutation_database.csv'

csv_file = args[1]
outPath = args[2]
sample_name = args[3]


df = pd.DataFrame()
df = pd.read_csv(csv_file)
df_mutDb = pd.read_csv(mutDb_file)

info = dict()
info.setdefault('Mutation', [])
info.setdefault('Genomic', [])
info.setdefault('Protein', [])
info.setdefault('Nucleotide', [])
info.setdefault('Comment', [])
info.setdefault('PMID', [])



for index, row in df.iterrows():

    match=False
    chrm=row['Chr']
    start=row['Start']
    end=row['End']
    ref=row['Ref']
    alt=row['Alt']
    
    for i,r in df_mutDb.iterrows():
        
        db_chrm=df_mutDb.loc[i,'Chr']
        db_start=df_mutDb.loc[i,'Start']
        db_end=df_mutDb.loc[i,'End']
        db_ref=df_mutDb.loc[i,'Ref']
        db_alt=df_mutDb.loc[i,'Alt']
        
        
        if (chrm==db_chrm) and (start==db_start) and (end==db_end) and (ref==db_ref) and (alt==db_alt):
            db_Mutation=df_mutDb.loc[i,'Mutation']
            db_info1=df_mutDb.loc[i,'Genomic']
            db_info2=df_mutDb.loc[i,'Protein']
            db_info3=df_mutDb.loc[i,'Nucleotide']
            db_info4=df_mutDb.loc[i,'Comment']
            db_info5=df_mutDb.loc[i,'PMID']         
            info['Mutation'].append(db_Mutation)
            info['Genomic'].append(db_info1)
            info['Protein'].append(db_info2)
            info['Nucleotide'].append(db_info3)
            info['Comment'].append(db_info4)
            info['PMID'].append(db_info5)
            break
        elif (i==len(df_mutDb)-1):
            info['Mutation'].append("")
            info['Genomic'].append("")
            info['Protein'].append("")
            info['Nucleotide'].append("")
            info['Comment'].append("")
            info['PMID'].append("")


df2=pd.DataFrame(info, columns=info.keys())
df_final = pd.concat([df, df2], axis=1)

df_final.to_csv(outPath+f'/{sample_name}_combined.csv', index=False)


