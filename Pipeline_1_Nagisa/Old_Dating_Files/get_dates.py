"""
Automated Genbank queries to grab collection dates when available
Requires lxml
"""
import os
import sys
from seqUtils import convert_fasta
import re
import gb
from time import sleep
from csv import DictWriter

# regex for matching accession numbers
p = re.compile('[bdegjm]+_([A-Z]+[0-9]+)|(NC_[0-9]+)|([A-Z]{1,2}[0-9]{5,6})[^0-9]|([A-Z]{1,2}[0-9]{5,6})$')

infile = open(sys.argv[1], 'rU')
fasta = convert_fasta(infile.readlines())
infile.close()

outfile = open(sys.argv[2], 'w')
writer = DictWriter(outfile, fieldnames=['accno', 'header', 'pmid', 'date', 'country', 'host', 
    'isolation_source', 'isolate', 'note', 'comment'], extrasaction='ignore')
writer.writeheader()

for i, (h, s) in enumerate(fasta):
    matches = p.findall(h)
    if len(matches) == 0:
        print ('unable to find accession number', h)
        continue
    
    acc = max(matches[0], key=len)
    print (acc)
    
    entry = gb.GenbankEntry(acc)
    if entry.xml is None:
        # failed to retrieve record
        outfile.write('%s,NA\n' % h)
        continue
    
    if not hasattr(entry, 'source_dict'):
        outfile.write('%s,,,,\n' % h)
        continue
    
    entry.source_dict.update({'accno': acc})
    entry.source_dict.update({'header': h})
    entry.source_dict.update({'comment': entry.comment})
    entry.source_dict.update({'pmid': '|'.join(entry.pmids)})
    
    writer.writerow(entry.source_dict)
    outfile.flush()
    
    sleep(0.5) # wait a second, to avoid spamming Genbank

outfile.close()

