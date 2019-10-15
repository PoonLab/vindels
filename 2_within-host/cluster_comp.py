import sys
import os
from subprocess import *
import pandas as pd  

if len(sys.argv) != 2:
	print("USAGE : python cluster_comp.py [desired node]")
	sys.exit()

target = sys.argv[1]

try:
	target = int(target)
	if target < 0 or target > 3: 
		print("Invalid node number")
		sys.exit()
except:	
	print("Invalid node number")
	sys.exit()

call = check_output(["bpstat", "-p"])
call = [x for x in call.split("\n") if x != '']

print(call)

df = pd.DataFrame({'process':[x.split('\t')[0] for x in call], 'node':[x.split('\t')[1] for x in call]})

print(df)

df = df.where(df['node'] == str(target))
df = df.dropna()
print(df)
