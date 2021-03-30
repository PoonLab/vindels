import os 
import sys
import pandas as pd
files = os.listdir(sys.argv[1])

pat = pd.Series([x.split("-")[0] for x in files])



print(pat.value_counts())
