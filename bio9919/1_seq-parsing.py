from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

data = list(SeqIO.parse("/home/jpalmer/Documents/bio9919-final-project/HCV-europe-1b-protF.fasta", "fasta"))

count = 0
countries = {}
with open("/home/jpalmer/Documents/bio9919-final-project/HCV-europe-1b-protF.fasta") as handle:
    for header,seq in SimpleFastaParser(handle):
        fields = header.split(".")
        
        if fields[0] == '-':
            count+=1
            #countries.update({fields[3]:header})
        count +=1
print(count)
print(countries)


#for record in data:
    #print(record.id, len(record))

#    header = record.id.split(".")
    #print(header)
