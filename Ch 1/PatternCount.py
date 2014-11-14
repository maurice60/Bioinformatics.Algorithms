__author__ = 'maurice'
from itertools import *

def fRead(fileName):
    txt = ""
    x = open(fileName, 'r')
    txt = x.read()
    x.close()
    return txt

def tRead(fileName):
    fd = open(fileName)
    st = ''.join([line.rstrip('\n').rstrip('\r') for line in fd if not line.startswith('>g')])
    fd.close()
    return st

def complementary(base):
    if  base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'

def findFreq(frequent):
    x = max(frequent.values())
    print 'Max freq', x
    out = []
    for j in frequent:
        if frequent[j] == x:
            out.append(j)
    return out

def patternCount(text, pattern):
    count = 0
    for i in range(0, len(text)-len(pattern)+1):
        if text[i:].startswith(pattern):
            count += 1
    return count

def patternWhere(text, pattern):
    idx = []
    for i in range(0, len(text)-len(pattern)+1):
        if text[i:].startswith(pattern):
            idx.append(i)
    return idx

def mostFrequent(text, k):
    frequent = {}
    for i in range(0, len(text) - k + 1):
        pat = text[i:i+k]
        if not frequent.has_key(pat):
            frequent[pat] = patternCount(text, pat)
    return findFreq(frequent)

def occurrances(text, k):
    frequent = {}
    for i in range(0, len(text) - k + 1):
        pat = text[i:i+k]
        if not frequent.has_key(pat):
            frequent[pat] = patternCount(text, pat)
    return frequent

def reverseComplement(text):
    ans = []
    for c in text:
        ans.append(complementary(c))
    ans.reverse()
    return ''.join(ans)

def clumpFind(text, k, l, t):
    out = []
    fDict = occurrances(text[0:l], k)
    for z in fDict:
        if fDict[z] >= t:
            out.append(z)
    for x in range(1, len(text) - l):
        fDict[text[x-1:x-1+k]] -= 1
        end = text[x:x+l][-k:]
        if fDict.has_key(end):
            fDict[end] += 1
            if fDict[end] >= t and end not in out:
                out.append(end)
        else:
            fDict[end] = 1
    return out

def skewCount(text):
    skew = 0
    ans = [0]
    for c in text:
        if c == 'G':
            skew += 1
        elif c == 'C':
            skew -= 1
        ans.append(skew)
    return ans

def minimumSkew(text):
    ans = []
    lis = skewCount(text)
    k = min(lis)
    for x, y in enumerate(lis):
        if y == k:
            ans.append(x)
    return ans

def maximumSkew(text):
    ans = []
    lis = skewCount(text)
    k = max(lis)
    for x, y in enumerate(lis):
        if y == k:
            ans.append(x)
    return ans

def hammingDistance(p, q):
    d = 0
    for i, c in enumerate(p):
        if q[i] != c:
            d += 1
    return d

def approxMatch(text, pattern, n):
    idx = []
    for i in range(0, len(text)-len(pattern)+1):
        if hammingDistance(text[i:i+len(pattern)], pattern) <= n:
            idx.append(i)
    return idx

def approximatePatternCount(text, pattern, n):
    count = 0
    for i in range(0, len(text)-len(pattern)+1):
        if hammingDistance(text[i:i+len(pattern)], pattern) <= n:
            count += 1
    return count

def allWords(k):
    out = set()
    base = ['A', 'C', 'G', 'T']
    ret = list(product(base, repeat = k))
    for t in ret:
        out.add(''.join(t))
    return out

def allPatterns(text, k):
    out = {}
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        out[pattern] = patternCount(text, pattern)
    return out

def hammingDistanceMax(p, q, d):
    h = 0
    for i, c in enumerate(p):
        if q[i] != c:
            h += 1
            if h > d:
                return False
    return True

def neighbours(pattern, d):
    if d == 0:
        return pattern
    return neighboursA(list(pattern), d)

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

def mostFrequentApproximate(text, k, d):
    frequent = {}
    for i in range(0, len(text) - k + 1):
        pat = text[i:i+k]
        for neigh in neighbours(pat, d):
            if frequent.has_key(neigh):
                frequent[neigh] += 1
            else:
                frequent[neigh] = 1
    return findFreq(frequent)

def mostFrequentApproximateReverse(text, k, d):
    frequent = {}
    for i in range(0, len(text) - k + 1):
        pat = text[i:i+k]
        for neigh in neighbours(pat, d):
            if frequent.has_key(neigh):
                frequent[neigh] += 1
            else:
                frequent[neigh] = 1
        pat = reverseComplement(pat)
        for neigh in neighbours(pat, d):
            if frequent.has_key(neigh):
                frequent[neigh] += 1
            else:
                frequent[neigh] = 1
    return findFreq(frequent)


# g = tRead('../Salmonella_enterica.txt')
# w = minimumSkew(g)
# print w
# for x in mostFrequentApproximateReverse(g[w[0]:w[-1]+500], 9, 1):
#     print(x),
# for x in mostFrequentApproximateReverse(g[w[0]:w[-1]+1000], 9, 1):
#     print(x),
# print(approximatePatternCount(fRead('dna.txt'), fRead('pat.txt'), 2))
# x = tRead('../Salmonella_enterica.txt')
# print skewCount('CATTCCAGTACTTCATGATGGCGTGAAGA')
# print len(neighbours('CCCC', 3))