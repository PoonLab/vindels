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

    for parent in root.iter():

        # removes all operators responsible for modifying the tree 
        if parent.tag == "operators":
            for toRemove in ["subtreeSlide", "narrowExchange", "wideExchange", "wilsonBalding"]:
                element = parent.find(toRemove)
                if element != None:
                    #print(element.tag)
                    parent.remove(element)

        # changes the population size to 5 
        if parent.tag == "populationSizes":
            if parent != None:
                element = parent.find("parameter")
                element.attrib = {"dimension":"5", "id":"skyline.popSize" ,"lower":"0.0" ,"value":"1.0"}

        if parent.tag == "groupSizes":
            if parent != None:
                element = parent.find("parameter")
                element.attrib = {"dimension":"5", "id":"skyline.groupSize"}


        # removes the simulator element and replaces it with the starting tree
        if parent.tag == "coalescentSimulator":
            parent.clear()
            parent.tag = "newick"
            parent.attrib = {'id':'startingTree'}
            parent.text = "\n"+tree
    
    xml.write("/home/jpalmer/PycharmProjects/hiv-withinhost/5_1_BEASTguided/"+xmlname)
    


