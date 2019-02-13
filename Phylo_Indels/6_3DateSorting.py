from glob import glob
from seqUtils import *
import calendar
import os
#uses the subtype dictionary.txt file to sort accession numbers + dates into the relevant output file

#dictionary for ABBREVIATED MONTH NAMES
caldict = {}
for k,v in enumerate(calendar.month_abbr):
    if len(str(k)) == 1:
        caldict[v] = "0"+str(k)
    else:
        caldict[v] = str(k)

#Dictionary for FULL MONTH NAMES
caldict2 = {}
for i, x in enumerate(calendar.month_name):
    if len(str(i)) == 1:
        caldict2[x] = "0"+str(i)
    else:
        caldict2[x] = str(i)
print(caldict2)


#---------------------------------------------

count = 0

input = open("genbank_dates.txt" , 'r')
dictionary = open("dictionary.txt", 'r')


gbdates = {}

#GENBANK DICTIONARY : [ACCNO] = genbank date
for line in input:
    info = line.strip("\n\r").split(".")

    gbdates[info[0]] = info[1]



# ACCDICT DICTIONARY :  [accnos] = (default year, subtype)
# SUBTYPES DICTIONARY :  [subtype] =  accession #s
# used to iterate through all accnos in one subtype file at a time
accdict = {}
subtypes = {}
for n in dictionary:
    data = n.strip("\n\r").split(",")
    accdict[data[0]] = data[1]
    if data[2] in subtypes.keys():
        subtypes[data[2]].append(data[0])

    else:
        subtypes[data[2]] = []
        subtypes[data[2]].append(data[0])

print()

for x in subtypes:
    print(x)

    #creates a temp file containing all the accurate genbank dates

    temp = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/Dates/temp.txt", 'w')
    for accno in subtypes[x]:    #where i is the accession number
        temp.write(accno + "," + gbdates[accno] + '\n')   #write the accession number and the genbank date associated with it
    temp.close()
    infile = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/Dates/temp.txt", 'r')

    output = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/Dates/" + x +".csv", 'w')


    #then read this date file; edit date formats and insert default dates if needed
    for n, line in enumerate(infile):

        info = line.strip("\n").split(",")
        accno = info[0]
        date = info[1]
        datefill = date

        if "/" in line:

            year = date.split("/")[1]

            date = year + "-01"
            print(line)
        # Handles the rare September 28, 2008 format
        elif len(info) == 3:

            dmy = date.split(" ")

            holder = dmy[1]

            dmy[1] = caldict2[dmy[0]]

            dmy[0] = info[2].strip(" ")
            if len(holder) == 2:
                dmy.append(holder)
            else:
                dmy.append("0" + holder)

            date = dmy[0] + "-" + dmy[1] + "-" + dmy[2]
            # datefill = date
        # Handles all normal dates with
        else:
            dmy = date.split("-")

            # 2008 format handling
            if len(dmy) == 1 and dmy[0].isdigit():

                date = dmy[0]

            # 'blank space' handling
            elif len(dmy) == 1 and dmy[0] == " ":
                year = accdict[accno]

                date = year

            # Aug-2008 format handling
            elif len(dmy) == 2:

                month = caldict[dmy[0]]

                del dmy[0]

                dmy.append(month)

                date = dmy[0] + "-" + dmy[1]

            # 13-Aug-2008 format handling
            elif len(dmy) == 3:

                holder = dmy[2]
                dmy[2] = dmy[0]
                dmy[0] = holder

                dmy[1] = caldict[dmy[1]]

                date = dmy[0] + "-" + dmy[1] + "-" + dmy[2]

        output.write(accno + "," + date + "\n")
    output.close()


    reopen = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/Dates/" + x + ".csv", 'r')
    newfile = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/Dates/edit/" + x + ".csv", 'w')

    for line in reopen:
        data = line.strip("\n").split(",")
        yr = data[1].split("-")[0]
        if int(yr) > 2017 or int(yr) < 1955:
            print(data[0])
            print(data[1])
        if data[0] == "EF593221":
            data[1] = "2001-11-30"
        elif data[0] == "EF593219":
            data[1] = "2006-08-31"

        newfile.write(data[0] + "," + data[1] + "\n")

os.remove("/home/jpalmer/PycharmProjects/hiv-evolution-master/Dates/temp.txt")