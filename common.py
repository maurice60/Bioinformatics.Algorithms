__author__ = 'maurice'
def fRead(fileName):
    with open(fileName, 'r') as x:
        txt = x.read()
    return txt

def tRead(fileName):
    fd = open(fileName)
    st = ''.join([line.rstrip('\n').rstrip('\r') for line in fd if not line.startswith('>g')])
    fd.close()
    return st

def lReadT(fileName):
    x = open(fileName, 'r')
    txt = x.read()
    x.close()
    out = []
    for i in txt.split(" "):
        out.append(i)
    return out

def hammingDistance(p, q):
    d = 0
    for i, c in enumerate(p):
        if q[i] != c:
            d += 1
    return d

def neighboursA(pattern, d):
    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']
    neighbourhood = set()
    suffixNeighbours = neighboursA(pattern[1:], d)
    for text in suffixNeighbours:
        if hammingDistance(pattern[1:], text)  < d:
            for x in ['A', 'C', 'G', 'T']:
                neighbourhood.add(x+text)
        else:
            neighbourhood.add(pattern[0]+text)
    return neighbourhood

def neighbours(pattern, d):
    if d == 0:
        return pattern
    return neighboursA(list(pattern), d)

