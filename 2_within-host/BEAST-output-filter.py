# script used to filter a BEAST output folder
# 1) removes any log files that are incomplete 
# 2) removes any other tree and operator files

import sys
import os
from glob import glob 
import shutil
import subprocess

inFolder = sys.argv[1]

if not inFolder.endswith("/"):
    inFolder += "/"

if not os.path.isdir(inFolder):
    print("this is not a proper directory")
    exit()
if not os.path.isdir(inFolder+"trees/"):
    os.mkdir(inFolder+"trees/")

if not os.path.isdir(inFolder+"operators/"):
    os.mkdir(inFolder+"operators/")

if not os.path.isdir(inFolder+"unfinished/"):
    os.mkdir(inFolder+"unfinished/")


files = glob(inFolder + "*")

for f in files: 


    if f.endswith(".log"):
        
        logFile = open(f, "rU")

        lastLine = subprocess.check_output(['tail', '-1', f])
        trace = lastLine.split()

        base = os.path.basename(f)
        tree = base.split(".")[0] + ".time.trees"
        ops = base.split(".")[0] + ".ops"

        if len(trace) > 0:

            if trace[0] != "100000000":
                os.rename(f, inFolder+"unfinished/"+base)
                if os.path.isfile(inFolder+tree):
                    os.rename(inFolder+tree, inFolder+"unfinished/"+tree)
                if os.path.isfile(inFolder+ops):
                    os.rename(inFolder+ops, inFolder+"unfinished/"+ops)
            
            else:
                if os.path.isfile(inFolder+tree):
                    os.rename(inFolder+tree, inFolder+"trees/"+tree)
                if os.path.isfile(inFolder+ops):
                    os.rename(inFolder+ops, inFolder+"operators/"+ops)

        else:
            os.rename(f, inFolder+"unfinished/"+base)
            if os.path.isfile(inFolder+tree):
                os.rename(inFolder+tree, inFolder+"unfinished/"+tree)
            if os.path.isfile(inFolder+ops):
                os.rename(inFolder+ops, inFolder+"unfinished/"+ops)


    
