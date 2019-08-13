# for checking whether study accession numbers are found within my data set 
from seqUtils import * 
import os 
import re 

acclist = [ "KC862610-KC863648", "HM638963-HM638965", "HM638988-HM638989", "HM638993-HM638994", "HM638997-HM639008", "HM639101-HM639103", "HM639105-HM639116", "HM639136-HM639139", "HM639143-HM639147", "HM639152-HM639157", "HM639159-HM639162", "HM639182-HM639184", "HM639199-HM639200", "HM639202-HM639203", "HM639078", "HM639084", "HM639096", "HM639141", "HM639150", "HM639168", "HM639171", "HM639173", "HM639177", "HM639180", "HM639192", "HM639197"]


acclist2 = [ 'EU574937-EU575065', 'EU575067-EU575212', 'EU575214-EU575231', 'EU575233','EU575235-EU575251', 'EU575253-EU575265',
'EU575267-EU575272', 'EU575274-EU575441', 'EU575443-EU575468', 'EU575470-EU575552', "EU575554-EU575636", "EU575638-EU575704", "EU575706-EU575775", "EU575777-EU575852", "EU575854-EU575943", "EU575945-EU575980", "EU575982-EU575990", "EU575992-EU576064", "EU576066-EU576089","EU576091-EU576237", "EU576239-EU576292", "EU576294-EU576296","EU576621-EU576642", "EU576644", "EU576646-EU576774", "EU576776-EU576799", "EU576801-EU576814",
"EU576816-EU576817", "EU576819-EU576840", "EU576842-EU576936", "EU576938-EU577005", "EU577007-EU577100", "EU577102-EU577114", 
"EU577116-EU577310", "EU577312-EU577350", "EU577352-EU577433", "EU577435-EU577440", "EU577442-EU577478", "EU577480-EU577662", 
"EU577664-EU578089", "EU578091-EU578109", "EU578111-EU578174", "EU578176-EU578239", "EU578241-EU578292", "EU578294-EU578307", 
'EU578309-EU578321', "EU578323-EU578328", "EU578330-EU578331", "EU578333-EU578375", "EU578377-EU578512", "EU578514-EU578559", 
"EU578561-EU578576", "EU578578", "EU578580-EU578636", "EU578638-EU578677", "EU578679-EU578686", "FJ495818", "FJ495937",
"FJ496000-FJ496001", "GU330247-GU330646", "GU330648-GU330861", "GU330938-GU331030", "GU331032-GU331095","GU331098-GU331102", "GU331114", "GU331116", "GU331121-GU331122",
"GU331127", "GU331129-GU331133", "GU331183-GU331217",
"GU331634", "GU331721"]

infile = open("/home/jpalmer/PycharmProjects/glyc-analysis/7_prnseq/gp120.fasta", "rU")

fasta = parse_fasta(infile)

allaccnos = [x.split(".")[1] for x in fasta]
#print(accnos[1])
total = 0
count = 0
for accno in acclist:
    if "-" in accno:
        accno1, accno2 = accno.split("-")
        start = int(accno1[2:])
        end = int(accno2[2:])
        prefix = accno1[:2]
        for x in range(start,end):
            total += 1
            accno = prefix + str(x)
            
            if accno in allaccnos:
                count+=1
                print(accno)
    
    else:
        total += 1
        if accno in allaccnos:
            count += 1
            print(accno)
print(total)
print(count)
    