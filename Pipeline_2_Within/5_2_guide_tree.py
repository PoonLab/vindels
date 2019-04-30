from Bio import Phylo
import cStringIO
from glob import glob 
import sys
import os
import xml.etree.ElementTree as ET 

xmlFolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/*.xml")

for infile in xmlFolder:
    xml = ET.parse(infile)
    xmlname = os.path.basename(infile)
    treename = "RAxML_bestTree."+ xmlname.split("-")[0] + ".tree"

    treefile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/4_5_Raxml/guide_trees/"+treename, "r")

    tree = treefile.readline()

    root = xml.getroot()

    for child in root:
        if child.tag == "coalescentSimulator":
            child.clear()
            child.tag = "newick"
            child.attrib = {'id':'startingTree'}

            print(child.tag, child.attrib)
            child.text = "\n"+tree

    xml.write("/home/jpalmer/PycharmProjects/hiv-withinhost/5_1_BEASTguided/"+xmlname)


