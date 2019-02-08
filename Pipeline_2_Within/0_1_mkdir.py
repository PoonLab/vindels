import re
from glob import glob
import os
from seqUtils import *



folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fastamsa")

for file in folder:
    filename = file.split("/")[-2]
    print(filename)
    #make directories
    #os.mkdir("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/" + filename)
    #rename
    newfile = file[:-3]
    print(newfile)
    os.rename(file, newfile)