from Bio import Phylo
import cStringIO
from glob import glob 
import sys
import os
import xml.etree.ElementTree as ET 
import lxml
import re 
'''
for s in range(len(sys.argv)):
    if not sys.argv[s].endswith("/"):
        sys.argv[s] = sys.argv[s] + "/"'''

xmlFolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/2e5-lognorm/*.xml")

#outpath = '/home/jpalmer/12BEAST/output/'

for infile in xmlFolder:
    xml = ET.parse(infile)
    xmlname = os.path.basename(infile) 
    treename = xmlname.split("-")[0] + ".tree"

    treefile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/4_5_Raxml/guide_trees/"+treename, "r")

    tree = treefile.readline()
    root = xml.getroot()

    seqcount = 0
    removelist = []
    #tempremove = []
    dates = []
    for element in root.iter():

        # POP SIZE = 5 
        # ---------------------
        # changes the population size to 5 
        '''if element.tag == "populationSizes":
            if element != None:
                elem = element.find("parameter")
                elem.set("dimension","5")
        if element.tag == "groupSizes":
            if element != None:
                elem = element.find("parameter")
                elem.set("dimension","5")'''

        # FIXING THE GUIDE TREE 
        # -----------------------------------------------
        # sets guide tree and removes all operators responsible for modifying the tree 
        if element.tag == "operators":
            for toRemove in ["subtreeSlide", "narrowExchange", "wideExchange", "wilsonBalding"]:
                elem = element.find(toRemove)
                if elem != None:
                    #print(elem.tag)
                    element.remove(elem)

        # sets the guide tree 
        if element.tag == "coalescentSimulator":
            element.clear()
            element.tag = "newick"
            element.attrib = {'id':'startingTree'}
            element.text = "\n"+tree
        
        
        '''# for editing the output folder path  
        # ------------------------------------------------       
        if element.get('id') == "mcmc":
            old = element.get('operatorAnalysis')
            element.set('operatorAnalysis', outpath+old.split("/")[-1])
            #print(element.get('operatorAnalysis'))
        '''
        
        # For Bayesian Skygrid Coalescent
        
        #if element.tag == "date":
        #    dates.append(float(element.get("value")))
        


        '''# For Bayesian skyride 
        if element.get("id") != None:
            tip_check = re.search('.*\..*\..*\..*\..*\..*\..*\_\d*', element.get("id"))
        else:
            tip_check = None
        if tip_check != None:
            seqcount += 1'''

        '''
        # STRICT CLOCK MODEL
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
        '''

    '''for r in removelist:
        #print(r)
        root.remove(r)'''

    # For editing the SkyGrid
    #for element in root.iter():
    #    if element.get("id") == "skygrid.cutOff":
    #        print(element.attrib)
    #        element.set("value", str(max(dates)))
    
    xml.write("/home/jpalmer/PycharmProjects/hiv-withinhost/5_1_BEASTguided/2e5-lognorm/"+xmlname)
    
    


