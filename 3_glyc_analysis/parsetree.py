import sys

from Bio import Phylo

tr = Phylo.read(sys.argv[1], 'newick')
target = int(sys.argv[2])

tips = tr.get_terminals()

if target >= len(tips):
  sys.exit()

while len(tips) > target:
  # find shortest tip
  lengths = [(tip.branch_length, tip) for tip in tips]
lengths.sort()
min_length, min_tip = lengths[0]
tr.prune(min_tip)
tips = tr.get_terminals()

Phylo.write(tr, sys.stdout, 'newick')