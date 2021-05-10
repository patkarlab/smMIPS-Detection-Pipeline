import pandas as pd
import sys

args = sys.argv

file1 = args[1]
file2 = args[2]
outfile = args[3]

df1 = pd.read_csv(file1, delimiter="\t")
df2= pd.read_csv(file2, delimiter="\t")

concat= pd.concat([df1,df2],axis=0)

concat.drop(columns=['ID'], inplace=True)
concat.replace(to_replace='.', value='-1', inplace=True)
concat.to_csv(outfile, index=False)
