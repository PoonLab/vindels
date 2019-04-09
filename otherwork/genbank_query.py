import sys
import re

def reformat_dates(dateList):
    outputDates = []
    for date in dateList:        
        
        # RARE FORMAT = 08/2008 
        '''if "/" in date:
            print(date)
            year = date.split("/")[1]

            date = year + "-01-01"
            print(year)

        #RARE FORMAT =  September 28, 2008
        elif " " in date:
            dateFields = date.split(" ")

            holder = dateFields[1]

            dateFields[1] = caldict2[dateFields[0]]

            dateFields[0] = info[2].strip(" ")
            if len(holder) == 2:
                dateFields.append(holder)
            else:
                dateFields.append("0"+holder)

            finalDate = dateFields[0] + "-" + dateFields[1] + "-" + dateFields[2]'''

        #Handles all normal dates with
        
        dateFields = date.split("-")

        #2008 format handling
        if len(dateFields) == 1 and dateFields[0].isdigit():
            finalDate = dateFields[0]


        #'blank space' handling
        '''elif len(dateFields) == 1 and dateFields[0] == " ":
            year = accdict[accno][0]
            finalDate = year'''


        #Aug-2008 format handling
        elif len(dateFields) == 2:
            finalDate = dateFields[1]


        #13-Aug-2008 format handling
        elif len(dateFields) == 3:
            finalDate = dateFields[2]
        outputDates.append(finalDate)
    return outputDates



gbFile = open("/home/jpalmer/Downloads/staph_aureus_resistance.gb",'r')


recordDictionary = {}
recordLines = []
for line in gbFile: 
    endOfRecord = re.search("\/\/\n", line)
    if endOfRecord != None:
        #print(line)
        #dump all the saved lines into a dictionary
    
        Found = False
        count=0

        #search for the accession number within each record
        while Found == False and count < len(recordLines):
            #print(count)
            if "ACCESSION" in recordLines[count]:
                Found = True
                accessionNo = recordLines[count].split()[1]
            count+=1
        if not Found:
            print("THIS SEQUENCE HAS A PROBLEM")
        #print(accessionNo)

        recordDictionary[accessionNo] = recordLines
        recordLines = []
    else:
        #otherwise add the current line to the saved set of lines
        recordLines.append(line)
#print(recordDictionary.keys())
#print(recordDictionary)

gbFile.close()

organismList = []
plasmidList = []
dateList = []
countryList = []
lenList = []
hostList = []
seqList = []

for record in recordDictionary:
    lenFound = False
    orgFound = False
    plasmidFound = False
    dateFound = False
    countryFound = False
    hostFound = False

    for line in recordDictionary[record]:   # reading through the lines of the record
        if "LOCUS" in line and "bp" in line:
            fields = line.split()
            seqLength = fields[2]
            lenList.append(seqLength)

        if "/organism=" in line:
            fields = line.split("=")
            organism = fields[1].strip('"\n')
            words = organism.split()
            if len(words) > 2:
                organism = " ".join(words[0:2])
            #print(organism)
            organismList.append(organism)
        if "/plasmid=" in line:
            fields = line.split("=")
            plasmid = fields[1].strip('"\n')
            if "unnamed" in plasmid or plasmid == "2":
                plasmid = "NA"
            plasmidList.append(plasmid)
        if "/collection_date=" in line:
            fields = line.split("=")
            colDate = fields[1].strip('"\n')
            
            dateList.append(colDate)
            print(colDate)
        
        #if line contains country, append country to countryList


        #if line contains host, append host to hostList



