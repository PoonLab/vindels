import re
from glob import glob
import os
from seqUtils import *

folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/conserved/*/")

for file in folder:
    filename = file.split("/")[-2]
    print(filename)
    os.mkdir("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/" + filename)