import sys
import pandas as pd

args = sys.argv
filename = args[1]
sample_name = args[2]
directory = args[3]
df = pd.read_csv(filename)

x = df['Otherinfo']


data = dict()

data.setdefault('REF_COUNT', [])
data.setdefault('ALT_COUNT', [])
data.setdefault('VAF', [])

for row in x:
    rowitems=row.split('\t')
    formatval=rowitems[-3].split(';')
    
    
    vaf=float(formatval[1].split('=')[1])  #AF
    dp4=formatval[3].split('=')[1]  #DP4
    ref=float(dp4.split(',')[0]) + float(dp4.split(',')[1])  #Ref forward+reverse
    alt=float(dp4.split(',')[2]) + float(dp4.split(',')[3])  #Alt forward+reverse
    
    data['REF_COUNT'].append(ref)
    data['ALT_COUNT'].append(alt)
    data['VAF'].append("{:.2%}".format(vaf))


df1=df.iloc[:,:5]
df2=pd.DataFrame(data, columns=data.keys())
df3=df.iloc[:,5:]

horizontal_stack = pd.concat([df1, df2, df3], axis=1)
horizontal_stack.rename(columns = {'Func.refGene':'Variant Site', 'ExonicFunc.refGene':'Variant Function'}, inplace = True)
outfile = directory + sample_name + "_lofreq" + ".csv"
horizontal_stack.to_csv(f'{outfile}', index=False)



