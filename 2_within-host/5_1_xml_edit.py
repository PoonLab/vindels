from Bio import Phylo
import cStringIO
from glob import glob 
import sys
import os
import xml.etree.ElementTree as ET 
import lxml

'''
for s in range(len(sys.argv)):
    if not sys.argv[s].endswith("/"):
        sys.argv[s] = sys.argv[s] + "/"'''

xmlFolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/strict/*.xml")

for infile in xmlFolder:
    xml = ET.parse(infile)
    xmlname = os.path.basename(infile)
    treename = "RAxML_bestTree."+ xmlname.split("-")[0] + ".tree"

    treefile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/4_5_Raxml/guide_trees/"+treename, "r")

    tree = treefile.readline()
    root = xml.getroot()
    previous = root

    removelist = []
    #tempremove = []

    for element in root.iter():
        #print(element.tag)


        # changes the population size to 5 
        if element.tag == "populationSizes":
            if element != None:
                elem = element.find("parameter")
                elem.set("dimension","5")

        if element.tag == "groupSizes":
            if element != None:
                elem = element.find("parameter")
                elem.set("dimension","5")

        # removes all operators responsible for modifying the tree 
        if element.tag == "operators":
            for toRemove in ["subtreeSlide", "narrowExchange", "wideExchange", "wilsonBalding"]:
                elem = element.find(toRemove)
                if elem != None:
                    #print(elem.tag)
                    element.remove(elem)
        if element.tag == "coalescentSimulator":
            element.clear()
            element.tag = "newick"
            element.attrib = {'id':'startingTree'}
            element.text = "\n"+tree
        


        # for changing the clock model to strict 
        # ------------------------------------------    
        if element.get("id") in ["coefficientOfVariation", "covariance"]:
            removelist.append(element)


        # remove the log elements 
        if element.get("id") == "fileLog":
            tempremove = []
            for elem in element.iter():
                if elem.get("idref") in ["coefficientOfVariation", "covariance", "ucld.stdev"]:
                    tempremove.append(elem)
            for i in tempremove:
                element.remove(i)

        if element.get("id") == "operators":
            tempremove = [element.find("swapOperator"), element.find("uniformIntegerOperator")]
            for e in element.iter():
                p = e.find("parameter")
                if p != None:
                    if p.get("idref") == "ucld.stdev":
                        tempremove.append(e)
            
            for i in tempremove:
                print(i.tag)
                element.remove(i)

        if element.get("id") == "prior":
            tempremove = []
            for e in element.iter():
                if e.get("mean") == "0.3333333333333333":
                    tempremove.append(e)
            
            for i in tempremove:
                element.remove(i)
    
        if element.get("idref") == "branchRates":
            element.tag = "strictClockBranchRates"

        if element.get("idref") == "ucld.mean":
            element.set("idref", "clock.rate")
        if element.get("label") == "ucld.mean":
            element.set("label", "clock.rate")

        if element.get("id") == "branchRates":
            element.clear()
            element.tag = "strictClockBranchRates"
            element.attrib = {'id':'branchRates'}
            element.append(ET.Element("rate"))
            element.find("rate").append(ET.Element("parameter", attrib={'id':'clock.rate', 'value':'1.0', 'lower':'0.0'}))

        previous = element
    for r in removelist:
        #print(r)
        root.remove(r)

    xml.write("/home/jpalmer/PycharmProjects/hiv-withinhost/5_2_strict/"+xmlname)
    


