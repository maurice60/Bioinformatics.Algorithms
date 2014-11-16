__author__ = 'maurice'
import sys
import os
# import math
# import random

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

# out = composition('0000111100101101000', 4, sort=True)
# out = genomePathWalk(aReadT('dna.txt'))
# out = overlapGraph(aReadT('dna.txt'))
# print out
# with open('outf.txt', 'w') as fil:
#     for x in out:
#         # print '{0} -> {1}'.format(*x)
#         fil.writelines('{0} -> {1}\n'.format(*x))
# a = 14
# l = []
# while True:
#     l.append(a)
#     if len(l) == 2**4:
#         break
#     a = (2*a) % (2**4)
#     while a in l:
#         a = (2*a + 1) % (2**4)
#
# bstr = ['{0:04b}'.format(d) for d in l]
#
# print genomePathWalk(bstr)
# print l, bstr
og = [(a,b) for (a,b) in overlapGraph(['{0:04b}'.format(d) for d in range(2**4)]) if a != b]



