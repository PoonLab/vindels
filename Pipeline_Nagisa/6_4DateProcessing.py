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
    for line in file:


        info = line.strip("\n").split(",")
        accno = info[0]
        date = info[1]
        datefill = date
        
        
        if "/" in line:
            print(line)
            year = date.split("/")[1]

            date = year + "-01-01"
            print(date)

        #Handles the rare September 28, 2008 format
        elif len(info) == 3:

            dmy = date.split(" ")

            holder = dmy[1]

            dmy[1] = caldict2[dmy[0]]

            dmy[0] = info[2].strip(" ")
            if len(holder) == 2:
                dmy.append(holder)
            else:
                dmy.append("0"+holder)

            date = dmy[0] + "-" + dmy[1] + "-" + dmy[2]
            #datefill = date
        #Handles all normal dates with
        else:
            dmy = date.split("-")

            #2008 format handling
            if len(dmy) == 1 and dmy[0].isdigit():


                date = dmy[0]
                #datefill = date + "-07-01"

            #'blank space' handling
            elif len(dmy) == 1 and dmy[0] == " ":
                year = accdict[accno][0]

                date = year
                #datefill = year + "-07-01"

            #Aug-2008 format handling
            elif len(dmy) == 2:

                month = caldict[dmy[0]]

                del dmy[0]


                dmy.append(month)

                date = dmy[0] + "-" + dmy[1]

                #datefill = date + "-15"

            #13-Aug-2008 format handling
            elif len(dmy) == 3:

                holder = dmy[2]
                dmy[2] = dmy[0]
                dmy[0] = holder

                dmy[1] = caldict[dmy[1]]

                date = dmy[0] + "-" + dmy[1] + "-" + dmy[2]

                #datefill = date
        output.write(accno + "," + date + "\n")
