# script used to filter a BEAST output folder
# 1) removes any log files that are incomplete 
# 2) removes any other tree and operator files

import sys
import os
from glob import glob 
from shutil import *
import subprocess

inFolder = sys.argv[1]

if not inFolder.endswith("/"):
    inFolder += "/"

if not os.path.isdir(inFolder):
    print("this is not a proper directory")
    exit()

runFolders = glob(inFolder+"*")

for folder in runFolders:

    finished = "/home/jpalmer/PycharmProjects/hiv-withinhost/8_1_Hfinished/"

    #if not os.path.isdir(unfinished):
        #os.mkdir(unfinished)

    files = glob(folder + "/*")

    if len(files) != 0:
        for f in files: 
            fileSize = subprocess.check_output(["stat", '--printf="%s"', f])
            fileSize = fileSize.strip('"')
            filename = os.path.basename(f)

            if int(fileSize) > 0:
                copyfile(f, finished+filename)


        
