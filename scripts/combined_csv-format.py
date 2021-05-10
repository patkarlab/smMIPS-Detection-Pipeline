import pandas as pd
import os, sys

args = sys.argv


varVcf = args[1]
mutectVcf = args[2]
lofreqVcf = args[3]
outPath = args[4]
sample_name = args[5]

files = [varVcf,mutectVcf,lofreqVcf]

df_var=pd.DataFrame()
df_mut=pd.DataFrame()
df_lofreq=pd.DataFrame()


df_var = pd.read_csv(varVcf)

df_mut = pd.read_csv(mutectVcf, sep=',')
df_lofreq = pd.read_csv(lofreqVcf, sep=',')
dataCaller = dict()
dataCaller.setdefault('Caller', [])

#ORDER DEPENDS ON ORDER OF FILES LIST
for i in range(len(df_var)):
    variantcaller=['VarDict']
    dataCaller['Caller'].append('VarDict')
for i in range(len(df_mut)):
    variantcaller=['Mutect2']
    dataCaller['Caller'].append('Mutect2')
for i in range(len(df_lofreq)):
    variantcaller=['Lofreq']
    dataCaller['Caller'].append('Lofreq')

df2 = pd.DataFrame()
for f in files:
    df2 = df2.append(pd.read_csv(f, sep=','), ignore_index=True)
df1=df2.iloc[:,:-1]
dfCaller=pd.DataFrame(dataCaller, columns=dataCaller.keys())

df = pd.concat([dfCaller, df1], axis=1)

num_of_callers=[]
rows_to_discard=[]
for index, row in df.iterrows():
    caller_counter=1
    caller=row['Caller']
    chrm=row['Chr']
    start=row['Start']
    end=row['End']
    ref=row['Ref']
    alt=row['Alt']
    #print("currently checking variant:",index,caller,chrm,start,end,ref,alt)
    for i in range(index+1,len(df)):
        this_caller=df.loc[i,'Caller']
        this_chrm=df.loc[i,'Chr']
        this_start=df.loc[i,'Start']
        this_end=df.loc[i,'End']
        this_ref=df.loc[i,'Ref']
        this_alt=df.loc[i,'Alt']
        #print(this_caller,this_chrm,this_start,this_end,this_ref,this_alt)
        if (chrm==this_chrm) and (start==this_start) and (end==this_end) and (ref==this_ref) and (alt==this_alt):
            caller=caller+','+this_caller
            df.at[index,'Caller']=caller
            rows_to_discard.append(i)
            caller_counter=caller_counter+1
    num_of_callers.append(caller_counter)

dataCaller_counter = dict()
dataCaller_counter.setdefault('Caller_Counter', num_of_callers)
dfCallerCounter=pd.DataFrame(dataCaller_counter, columns=dataCaller_counter.keys())

df_final = pd.concat([dfCallerCounter, df], axis=1)

df_final.drop(rows_to_discard,inplace=True)
df_final.to_csv(f'{outPath}/{sample_name}_combined.csv', index=False)








