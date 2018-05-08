import sys


#processes the two separate sequence files retrieved from Genbank
#reads line by line, finds and prints accession number, finds and prints collection date if it exists


# these files contain the raw data from GenBank
file1 = open("sequence1.txt", 'r')
file2 = open("sequence2.txt", 'r')

count = 0

dict = {}

date = ''

hasdate = True

#prints the date to this
output = open("genbank_dates.txt", 'w')

for line in file1:

    if "ACCESSION" in line:

        # hasdate will be false at this point if no date was located for the previous sequence
        # used to skip to the next line for the next sequence
        if hasdate == False:
            output.write(" \n")


        hasdate = False

        info = line.split(" ")
        accno = info[3].strip("\n") + "."
        output.write(accno)



    elif "/collection_date=" in line:
        data = line.split('=')
        date = data[1].strip('"\n') + "\n"

        hasdate = True
        output.write(date)


hasdate = True


for line in file2:

    if "ACCESSION" in line:

        if hasdate == False:
            output.write(" \n")

        hasdate = False

        info = line.split(" ")
        accno = info[3].strip("\n") + "."
        output.write(accno)



    elif "/collection_date=" in line:
        data = line.split('=')
        date = data[1].strip('"\n') + "\n"
        hasdate = True
        output.write(date)

output.close()