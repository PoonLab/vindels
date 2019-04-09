import sys
from glob import glob
from seqUtils import *

import calendar

#dictionary for ABBREVIATED MONTH NAMES
caldict = {}
for k,v in enumerate(calendar.month_abbr):
    if len(str(k)) == 1:
        caldict[v] = "0"+str(k)
    else:
        caldict[v] = str(k)
print(caldict)
#Dictionary for FULL MONTH NAMES
caldict2 = {}
for i, x in enumerate(calendar.month_name):
    if len(str(i)) == 1:
        caldict2[x] = "0"+str(i)
    else:
        caldict2[x] = str(i)
print(caldict2)

dictionary = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/dictionary.txt")
accdict = {}

#key = Accession Number ; Value = (Year, Subtype)
for line in dictionary:
    data = line.strip("\n\r").split(",")
    accdict[data[0]] = (data[1],data[2])


folder = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/Dates/*.txt")
for path in folder:

    file = open(path, 'r')

    #create the corresponding output file
    name = path.split("/")[-1].split(".")[0]
    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/Dates_edit/" + name + "2.txt",'w')


    #read through each line and edit the dates
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
            '''month = caldict[dateFields[0]]
            del dateFields[0]
            dateFields.append(month)'''
            finalDate = dateFields[1]


        #13-Aug-2008 format handling
        elif len(dateFields) == 3:
            '''holder = dateFields[2]
            dateFields[2] = dateFields[0]
            dateFields[0] = holder
            dateFields[1] = caldict[dateFields[1]]'''
            finalDate = dateFields[2]
        outputDates.append(finalDate)
    return outputDates
