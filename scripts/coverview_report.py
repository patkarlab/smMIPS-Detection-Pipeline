import pandas as pd
import os, sys

args = sys.argv

path = args[1]
outfile = args[2]
files = os.listdir(path)

dataframes=[]
commoncols=[]

for files in files:    
    path1, filename = os.path.split(files)
    samplename=os.path.splitext(filename)[0].split('.')[0]
    csvFile = f'{path}/{files}'
    df= pd.read_csv(csvFile, index_col=0)
    df=df.iloc[:,:6]
    df.rename(columns = {'Read count':samplename+'_Read count', 'Median coverage':samplename+'_Median coverage', 'Pass_or_flag':samplename+'_FILTER'},inplace=True)
    dataframes.append(df)
    
commoncols = list(dataframes[0].columns)[:3]
for i in range(1,len(dataframes)):
    dataframes[i].drop(commoncols, axis=1, inplace=True)
        
concat=pd.concat(dataframes, axis=1)
concat.to_csv(f'{outfile}/Coverview-Report.csv', index=None)
