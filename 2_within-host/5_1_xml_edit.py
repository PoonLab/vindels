from Bio import Phylo
import cStringIO
from glob import glob 
import sys
import os
import xml.etree.ElementTree as ET 
import lxml
import re 
import pandas as pd

for s in range(len(sys.argv)):
    if not sys.argv[s].endswith("/"):
        sys.argv[s] = sys.argv[s] + "/"

if len(sys.argv) != 2:
    print("USAGE: python 5_1_xml_edit.py [unique run id]")
    sys.exit()

run_id = sys.argv[1]

# read in the separate CSV containing all the root heights of the trees
df = pd.read_csv("/home/jpalmer/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/root-heights.csv", index_col="file")
print(df)

if not os.path.isdir('/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/' + run_id + "/"):
    print("ERROR: Preliminary run folder not found.")
    sys.exit()

if os.path.isdir('/home/jpalmer/PycharmProjects/hiv-withinhost/5_1_BEASTguided/' + run_id + "/"):
    print("ERROR: Run already exists as folder in 5_1. Please check the folder.")
    sys.exit()
else:
    os.mkdir('/home/jpalmer/PycharmProjects/hiv-withinhost/5_1_BEASTguided/' + run_id + "/")

# iterate through the 5BEAST/ --- folder directory    
xmlFolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/"+run_id+"/*.xml")

for infile in xmlFolder:
    xml = ET.parse(infile)
    xmlname = os.path.basename(infile) 

    if re.search("-original",xmlname) != None:
        outname = re.sub("-original","",xmlname)
    else:
        outname = xmlname
    
    print(xmlname)
    
    treename = xmlname.split("-")

    if len(treename) == 2: 
        treename = treename[0] + ".tree"
    else:
        treename = treename[0] + "-" +treename[1] + ".tree"

    # open the tree file 
    treefile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/guide_trees/"+treename, "rU")

    tree = treefile.readline()
    root = xml.getroot()

    seqcount = 0
    removelist = []
    #tempremove = []
    dates = []
    
    for element in root.iter():
        #print(element.tag)
        

        break    # used to avoid fixing the tree for BEAST run 53
        #if element.tag == "date":
        #   dates.append(float(element.get("value")))
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

        '''# SAMPLING FROM PRIOR ONLY 
        # -----------------------------------------
        if element.tag == "alignment":

            #<taxon idref="C.ZM.2002.16362.GQ485317.12.3_11"/>
            for seq in element.findall("sequence"):
                name = seq.find("taxon").get("idref")
                seq.clear()
                seq.append(ET.Element("taxon",{"idref":name}))
                seq.text="???"'''
                

        
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
        


        '''# STRICT CLOCK MODEL
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
        

    for r in removelist:
        #print(r)
        root.remove(r)'''

    # extracts the root_height value from the root height CSV file 
    rheight = df.loc[treename.split(".")[0],'root_height']
    print(rheight)

    # searches for the skygrid cutoff value and replaces it with the extracted root height
    for element in root.iter():
        if element.get("id") == "skygrid.cutOff":
            #print(element.attrib)
            element.set("value", str(rheight))
            print(element.attrib)
    xml.write("/home/jpalmer/PycharmProjects/hiv-withinhost/5_1_BEASTguided/"+run_id+"/"+outname)    



'''# for editing the output folder path  
# ------------------------------------------------       
if element.get('id') == "mcmc":
    old = element.get('operatorAnalysis')
    element.set('operatorAnalysis', outpath+old.split("/")[-1])
    #print(element.get('operatorAnalysis'))
'''




    


