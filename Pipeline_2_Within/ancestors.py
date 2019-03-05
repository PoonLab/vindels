from glob import * 
import sys 
import os 
from seqUtils import * 

def opposite(num):
    if num == 0:
        return 1
    else:
        return 0

def extractAncestors(inputFile):

    with open(inputFile) as handle:
        data = parse_fasta(handle)
    
    cherry = re.compile("\([^()]*:[^()]*,[^()]*\)")
    
    nodes = {}
    
    for header in data.keys():
        #print(header)
        check = cherry.match(header)

        if check != None:
            s1, s2  = header.strip("()").split(",")
            s1 = s1.split(":")[0]
            s2 = s2.split(":")[0]
            label =  s1 + "." + s2
            
            nodes[label] = [data[s1],data[s2],data[header]]
    

    iDict = {}
    dDict = {}
    for header in nodes.keys():
        #print(header)
        head1, head2 = header.split('.')
        for i in range(2):
            iTemp = ''
            dTemp = ''
            insertions = []
            deletions = []

            for n, Achar in enumerate(nodes[header][2]):
            
                Schar = nodes[header][i][n]
                Ochar = nodes[header][opposite(i)][n]

                #1 : ancestral gap, sequence gap  = nothing 
                #2 : ancestral gap, sequence char = insertion 
                #3 : ancestral char, sequence gap  = deletion 
                #4 : ancestral char, sequence char  = nothing 
                
                if Achar == "-":
                    #insertion in seq 1 
                    if Schar != '-':
                        iTemp += Schar
                        #print(Achar + "" + Schar+"" + Ochar)
                    #clear the dTemp 
                        if dTemp:
                            deletions.append(tuple((dTemp, n-1)))
                            dTemp = ''
                
                        
                    #nothing -- both gaps 
                    else:
                        #clear iTemp and dTemp
                        if iTemp:
                            insertions.append(tuple((iTemp, n-1)))
                            
                            iTemp = ''
                        if dTemp:
                            deletions.append(tuple((dTemp, n-1)))
                            dTemp = ''
                        
                elif Achar == "*":
                    #deletion in seq 1
                    if Schar == "-":
                        dTemp += Ochar
                        #print(Achar + "" + Schar+"" + Ochar)
                        #clear iTemp
                        if iTemp:
                            insertions.append(tuple((iTemp, n-1)))
                            iTemp = ''
                        
                    #nothing -- both have character
                    else:
                        #clear iTemp and dTemp
                        if iTemp:
                            insertions.append(tuple((iTemp, n-1)))
                            iTemp = ''
                        if dTemp:
                            deletions.append(tuple((dTemp, n-1)))
                            dTemp = ''
            print(insertions)
            print(deletions)
            iDict[] = insertions
            iDict[

        

def main():
    files = glob('/home/jpalmer/historian/bin/*.fasta')

    for f in files:
        extractAncestors(f)

if __name__ == '__main__': 
    main()

