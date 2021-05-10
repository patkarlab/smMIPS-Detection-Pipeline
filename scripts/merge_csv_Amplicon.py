import pandas as pd
import os
import sys

args = sys.argv
sample_name = args[1]
directory = args[2]

tsv_file=f'{directory}/{sample_name}/coverage/{sample_name}.counts.bed'
csv_table=pd.read_table(tsv_file,sep='\t',names=["Start", "End","Sequence","Reads"])
csv_table.to_csv(f'{directory}/{sample_name}/VariantCalling/{sample_name}coverage.counts.csv',index=False)

csvfile=[directory+f'/{sample_name}/VariantCalling/{sample_name}_combined.csv', directory+f'/{sample_name}/VariantCalling/{sample_name}coverage.counts.csv']
writer = pd.ExcelWriter(directory+f'/Final_Output/{sample_name}/{sample_name}_Final.xlsx')

for filename in csvfile:
    sheetname=os.path.split(filename)[1]
    df = pd.read_csv(filename)
    #print('process file:', filename, 'shape:', df.shape)
    df.to_excel(writer,sheet_name=os.path.splitext(sheetname)[0], index=False)
writer.save()
