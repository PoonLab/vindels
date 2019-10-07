from seqUtils import * 
import os 
import sys
import time
import glob


# this script will remove columns from an alignment that do not contain any character spaces 

# read in a file from standard input 

# convert fasta

# transpose fasta 

# iterate through the transposed fasta columns and check if any(chars) != "-"

# only add these columns to a new transposed fasta as output

# print to standard output 

def removeGaps(input_fasta):
    with open(input_fasta, "rU") as handle:
        tpfasta = transpose_fasta(convert_fasta(handle))

    filename = os.path.basename(input_fasta).split(".")[0]
    filedir = os.path.dirname(input_fasta)
    # generate a whitelist of all 
    whitelist = []
    for n, column in enumerate(tpfasta):
        
        #gaps = column.count("-")
        #freq = float(gaps)/len(x)
        #if freq < 0.95:
        #    whitelist.append(pos)
        if not all([x == "-" for x in column]):
            whitelist.append(n)
    
    with open(sys.argv[1], "rU") as handle:
        data = parse_fasta(handle)

    output = open(filedir+"/"+filename+"_nogaps.fasta", "w+")
    for header in data:
        newseq = ''
        
        for n, char in enumerate(data[header][:]):
            if n in whitelist:
                newseq += char
            
        output.write(">"+header+"\n"+newseq+"\n")
    output.close()
def main():
    if len(sys.argv) != 2:
        print("USAGE: python removeGaps.py [input FASTA file]")
        sys.exit()
    
    removeGaps(sys.argv[1])

        
        
        




if __name__ == '__main__': 
    main()
