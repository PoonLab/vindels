import sys


#processes the two separate sequence files retrieved from Genbank
#reads line by line, finds and prints accession number, finds and prints collection date if it exists


# these files contain the raw data from GenBank

dict = {}

date = ''



#prints the date to this
output = open("genbank_dates.txt", 'w')
for i in range(2):
    file = open('/home/jpalme56/PycharmProjects/hiv-evolution-master/gbdata'+str(i+1)+'.txt', 'r')

    hasdate = True
    for line in file:

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

output.close()