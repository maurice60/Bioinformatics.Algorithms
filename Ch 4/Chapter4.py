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

def deBruijnPrint(patterns):

    def deBruijn(patterns):
        out = sorted([(i[:-1], i[1:]) for i in patterns])
        return out

    lines = []
    last = 0
    lastUsed = ''
    for x in deBruijn(patterns):
        if x[0] == lastUsed:
            lines[last-1] += ',{0}'.format(x[1])
        else:
            lines.append('{0} -> {1}'.format(*x))
            last += 1
            lastUsed = x[0]
    return lines

# a = 15
# l = []
# while True:
#     l.append(a)
#     if len(l) == 2**4:
#         break
#     asf = a
#     a = (2*a) % (2**4)
#     while a in l:
#         a = (2*asf + 1) % (2**4)
#         asf += 1
#
bstr = ['{0:04b}'.format(d) for d in xrange(16)]
#
# print genomePathWalk(bstr)
# print l, bstr
# og = [(a,b) for (a,b) in overlapGraph(['{0:04b}'.format(d) for d in range(2**4)]) if a != b]
# print og
#
# out = composition('0000111100101101000', 4, sort=True)
# out = genomePathWalk(aReadT('dna.txt'))
# out = overlapGraph(aReadT('dna.txt'))
out = deBruijnPrint(bstr) #composition(fRead('dna.txt'), 12))
# print out
for x in out:
    print x
with open('outf.txt', 'w') as fil:
    fil.writelines('\n'.join(out))

