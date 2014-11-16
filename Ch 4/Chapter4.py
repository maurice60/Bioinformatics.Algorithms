__author__ = 'maurice'
import sys
import os
import math
import random

sys.path.append(os.path.relpath('../'))
from common import *

def genomePathWalk(dna):
    genome = dna[0]
    for x in dna[1:]:
        genome += x[-1]
    return genome

def overlapGraph(dna):
    out = []
    for s in dna:
        for h in [x for x in dna if x.startswith(s[1:])]:
            out.append((s, h))
    return sorted(out)

# out = composition(fRead('dna.txt'), 100, sort=True)
# out = genomePathWalk(aReadT('dna.txt'))
out = overlapGraph(aReadT('dna.txt'))
# print out
with open('outf.txt', 'w') as fil:
    for x in out:
        # print '{0} -> {1}'.format(*x)
        fil.writelines('{0} -> {1}\n'.format(*x))