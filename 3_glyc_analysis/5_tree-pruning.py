from seqUtils import *
from subprocess import Popen, PIPE

infile = open("/home/jpalmer/PycharmProjects/glyc-analysis/5_tree/RAxML_bestTree.gp120.tree", "rU")


p = Popen(" ".join(['python', '/home/jpalmer/vindels/3_glyc_analysis/prunetree.py', "/home/jpalmer/PycharmProjects/glyc-analysis/5_tree/RAxML_bestTree.gp120.tree", ">", "/home/jpalmer/PycharmProjects/glyc-analysis/6_pruned/pruned.tree" ]),shell=True)
