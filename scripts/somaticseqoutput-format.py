import pandas as pd
import sys

args = sys.argv

filename = args[1]
outfile = args[2]
vartools=['MuTect2, ', 'VarScan2, ', 'VarDict, ', 'LoFreq, ', 'Strelka']

df = pd.read_csv(filename)
x = df['Otherinfo1']
discarded_column=df.columns.get_loc('Otherinfo2')
data = dict()
somatic_cols=['Chr','Start','End','Ref','Alt','Variant_Callers','FILTER','SOMATIC_FLAG','VariantCaller_Count','REF_COUNT','ALT_COUNT','VAF','Func.refGene','Gene.refGene','ExonicFunc.refGene','AAChange.refGene','Gene_full_name.refGene','Function_description.refGene','Disease_description.refGene','cosmic84','PopFreqMax','1000G_ALL','ExAC_ALL','CG46','ESP6500siv2_ALL','InterVar_automated']
data.setdefault('FILTER', [])
data.setdefault('SOMATIC_FLAG', [])
data.setdefault('VariantCaller_Count', [])
data.setdefault('Variant_Callers', [])
data.setdefault('VAF', [])
data.setdefault('REF_COUNT', [])
data.setdefault('ALT_COUNT', [])
for row in x:
    rowitems=row.split('\t')
    data['FILTER'].append(rowitems[9])
    info=rowitems[10].split(';')
    formatval=rowitems[-1].split(':')
    readcounts=formatval[2]
    if len(info)!=5 and info[0]!='SOMATIC':
        data['SOMATIC_FLAG'].append('NON SOMATIC')
        name_toolskey=0
    else:
        data['SOMATIC_FLAG'].append(info[0])
        name_toolskey=1
    num_toolskey=name_toolskey+1    
    tools_binary=info[name_toolskey].split('=')[1].split(',')
    tools_name=''
    for key in range(len(tools_binary)):
        if tools_binary[key]=='1':
            tools_name+=vartools[key]
    data['Variant_Callers'].append(tools_name)        
    data['VariantCaller_Count'].append(info[num_toolskey].split('=')[1])
    vaf="{:.2%}".format(float(info[-1].split('=')[1]))
    data['VAF'].append(vaf)
    data['REF_COUNT'].append(readcounts.split(',')[1])
    data['ALT_COUNT'].append(readcounts.split(',')[-1])
            
df1=df.iloc[:,:5]
df2=pd.DataFrame(data, columns=data.keys())
df3=df.iloc[:,5:discarded_column]

horizontal_stack = pd.concat([df1, df2, df3], axis=1)
horizontal_stack.to_csv(outfile, index=False)
