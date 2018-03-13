import sys
from glob import glob
from seqUtils import *

import calendar

#dictionary for ABBREVIATED MONTH NAMES
dict = {}
for k,v in enumerate(calendar.month_abbr):
    if len(str(k)) == 1:
        dict[v] = "0"+str(k)
    else:
        dict[v] = str(k)


#Dictionary for FULL MONTH NAMES
dict2 = {}
for i, x in enumerate(calendar.month_name):
    if len(str(i)) == 1:
        dict2[x] = "0"+str(i)
    else:
        dict2[x] = str(i)


folder = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/Dates/*.txt")

for path in folder:

    file = open(path, 'r')

    #create the corresponding output file
    name = path.split("Dates")
    output = open(name[0] + "Dates_edit" + name[1],'w')

    
    #read through each line and edit the dates
    for line in file:
        info = line.strip("\n").split(",")
        accno = info[0]
        date = info[1]


        #Handles the rare September 28, 2008 format
        if len(info) == 3:
            print(info)
            dmy = date.split(" ")

            holder = dmy[1]

            dmy[1] = dict2[dmy[0]]

            dmy[0] = info[2].strip(" ")
            dmy.append(holder)

            date = dmy[0] + "-" + dmy[1] + "-" + dmy[2]

            output.write(accno + "," + date + "\n")

        #Handles all normal dates with
        else:
            dmy = date.split("-")

            #2008 format handling
            if len(dmy) == 1 and dmy[0].isdigit():
                dmy.append('07')
                dmy.append('01')

                date = dmy[0] + "-" + dmy[1] + "-" + dmy[2]

            #'blank space' handling
            elif len(dmy) == 1 and dmy[0] == " ":

                date = " "

            #Aug-2008 format handling
            elif len(dmy) == 2:

                month = dict[dmy[0]]

                del dmy[0]

                #Defaulted missing day to 15
                dmy.append(month)
                dmy.append('15')

                date = dmy[0] + "-" + dmy[1] + "-" + dmy[2]

            #13-Aug-2008 format handling
            elif len(dmy) == 3:

                holder = dmy[2]
                dmy[2] = dmy[0]
                dmy[0] = holder

                dmy[1] = dict[dmy[1]]

                date = dmy[0] + "-" + dmy[1] + "-" + dmy[2]

            output.write(accno + "," + date + "\n")
