import pandas as pd
import sys

args = sys.argv

file = args[1]
outfile = args[2]
regions = pd.read_csv(file, delimiter="\t")
					
regions.rename(columns = {'RC':'Read count', 'MEDCOV':'Median coverage', 'MINCOV':'Minimum coverage',
                          'MEDQCOV':'Median quality coverage','MINQCOV':'Minimum quality coverage',
                          'MAXFLMQ':'Maximum fraction of low mapping quality',
                          'MAXFLBQ':'Maximum fraction of low base quality'}, inplace = True) 
regions.to_csv(outfile, index=False)