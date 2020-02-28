import re 


infile = open("/home/jpalmer/slip-model.csv", "rU")
outfile = open("/home/jpalmer/slip-model-parsed.csv" , "w+")
newlines = []
col = []
for n, line in enumerate(infile):
    if n % 3 == 2:
        col.append(line.split(" ")[1].strip("\n"))

    if re.search("\[\d\] -\d+\n", line) == None:
        newline = re.sub('\[\d\] "STATE ', "", line)
        newline = re.sub(": ", "", newline)
        newline = re.sub('"', "", newline)
        newline = re.sub(" ", ",", newline)
        newlines.append(newline.strip("\n"))
print(newlines)
print(len(newlines))

for i in range(len(newlines)-1):
    print(i)
    outfile.write(newlines[i]+","+col[i]+"\n")