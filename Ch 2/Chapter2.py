__author__ = 'maurice'
import sqlite3
from collections import Counter
from operator import itemgetter

def reverseComplement(text):
    ans = []
    cx = sqlite3.connect("../Bioinformatics.db")
    px = cx.cursor()
    for c in text:
        px.execute("SELECT complement FROM complement WHERE base = ?", (c,))
        ans.append(px.fetchone()[0].encode('ascii','ignore'))
    cx.close()
    ans.reverse()
    return ''.join(ans)

def trimerFind(tri):
    cx = sqlite3.connect("../Bioinformatics.db")
    px = cx.cursor()
    px.execute("SELECT protein FROM codon WHERE trimer = ?", (tri,))
    out = px.fetchone()
    cx.close()
    if out is None:
        return ''
    return out[0]

def weightFind(base):
    cx = sqlite3.connect("../Bioinformatics.db")
    px = cx.cursor()
    px.execute("SELECT mass FROM moleculeMass WHERE molecule = ?", (base,))
    out = px.fetchone()
    cx.close()
    if out is None:
        return ''
    return out[0]

def allWeight():
    cx = sqlite3.connect("../Bioinformatics.db")
    px = cx.cursor()
    px.execute("SELECT mass FROM moleculeMass")
    out = px.fetchall()
    cx.close()
    ans = set()
    for x in out:
        ans.add(x[0])
    return ans

def aminoFind(ami):
    cx = sqlite3.connect("../Bioinformatics.db")
    px = cx.cursor()
    px.execute("SELECT trimer FROM codon WHERE protein = ?", (ami,))
    out = px.fetchall()
    cx.close()
    ans = []
    for x in out:
        ans.append(x[0].encode('ascii','ignore'))
    return ans

def fRead(fileName):
    x = open(fileName, 'r')
    txt = x.read()
    x.close()
    return txt

def tRead(fileName):
    fd = open(fileName)
    st = ''.join([line.rstrip('\n').rstrip('\r') for line in fd if not line.startswith('>g')])
    fd.close()
    return st

def lRead(fileName):
    txt = ""
    x = open(fileName, 'r')
    txt = x.read()
    x.close()
    out = []
    for i in txt.split(" "):
        out.append(int(i))
    return out

def lReadT(fileName):
    txt = ""
    x = open(fileName, 'r')
    txt = x.read()
    x.close()
    out = []
    for i in txt.split(" "):
        out.append(i)
    return out

def translation(text):
    out = ''
    for i in range(0, len(text), 3):
        out += trimerFind(''.join(text[i:i+3]))
    return out

def peptideNumber(peptide):
    ans = 1
    for c in peptide:
        n = len(aminoFind(c))
        print c, n
        ans *= n
    return ans

def peptideEncoding(text, peptide):
    l = []
    lr = []
    ans = []
    for c in peptide:
        a = aminoFind(c)
        ar = []
        ra = []
        for x in a:
            y = x.replace('U', 'T')
            ar.append(y)
            ra.append(reverseComplement(y))
        l.append(ar)
        lr.append(ra)
    lr.reverse()
    for i in range(len(text) - 3):
        found = True
        for j in range(len(l)):
            k = i + 3*j
            if text[k:k+3] not in l[j]:
                found = False
                break
        if found:
            ans.append(text[i:i+len(l)*3])
        found = True
        for j in range(len(lr)):
            k = i + 3*j
            if text[k:k+3] not in lr[j]:
                found = False
                break
        if found:
            ans.append(text[i:i+len(l)*3])
    return ans

def theoreticalSpectrum(text):
    ans = [0]
    t = 0
    for i, c in enumerate(text):
        x = weightFind(c)
        t += x
        ans.append(x)
        for j in range(len(text) - 2):
            x += weightFind(text[(i+j+1)%len(text)])
            ans.append(x)
    ans.append(t)
    return sorted(ans)

def linearSpectrum(text):
    ans = [0]
    t = 0
    for i, c in enumerate(text):
        x = weightFind(c)
        t += x
        ans.append(x)
        for j in range(1, len(text) - i):
            x += weightFind(text[i+j])
            ans.append(x)
    return sorted(ans)

def cyclopeptideSpectrum(lis):
    ans = [0]
    t = 0
    for i, x in enumerate(lis):
        t += x
        ans.append(x)
        for j in range(len(lis) - 2):
            x += lis[(i+j+1)%len(lis)]
            ans.append(x)
    ans.append(t)
    return sorted(ans)

def linearListSpectrum(lis):
    ans = [0]
    t = 0
    for i, x in enumerate(lis):
        t += x
        ans.append(x)
        for j in range(1, len(lis) - i):
            x += lis[i+j]
            ans.append(x)
    return sorted(ans)

def linearNumber(n):
    ans = 1
    for i in range(n, 0, -1):
        ans += i
    return ans

def peptideMass(p):
    m = 0
    for x in p:
        m += x
    return m

def consistent(spectrum, test):
    for t in test:
        if t not in spectrum:
            return False
        spectrum.remove(t)
    return True

def peptideScore(sp, test):
    def testDict(te):
        tp = {}
        for t in test:
            if tp.has_key(t):
                tp[t] += 1
            else:
                tp[t] = 1
        return tp
    x = 0
    tp = testDict(test)
    for t in tp.keys():
        if sp.has_key(t):
            x += min(sp[t], tp[t])
    return x

def cyclopeptideSequencing(spectrum):
    ans = []
    targ = max(spectrum)
    peptides = [[]]
    pLis = allWeight()
    while len(peptides) > 0:
        new = []
        for p in peptides:
            for l in pLis:
                x = p + [l]
                new.append(x)
        peptides = []
        for p in new:
            if peptideMass(p) == targ:
                if cyclopeptideSpectrum(p) == spectrum:
                    ans.append(p)
            if consistent(spectrum[:], p):
                peptides.append(p)
    return ans

def trim(toTrim, spectrum, n):
    scores = []
    for j in toTrim:
        scores.append((j, peptideScore(spectrum, linearListSpectrum(j))))
    srt = sorted(scores, key=itemgetter(1), reverse=True)
    if len(srt) >= n:
        targ = srt[n-1][1]
    else:
        targ = 0
    out = []
    for v in srt:
        if v[1] >= targ:
            out.append(v[0])
        else:
            break
    return out

def trimT(leaderboard, spectrum, n):
    toTrim = []
    pLis = []
    out = []
    for peptide in leaderboard:
        p = []
        for c in peptide:
            p.append(weightFind(c))
        toTrim.append(p)
        pLis.append((peptide, p))
    trimmed = trim(toTrim, Counter(spectrum), n)
    for lis in trimmed:
        for x in pLis:
            if x[1] == lis:
                out.append(x[0])
                break
    return out

def leaderboardCyclopeptideSequencing(spectrum, n, pLis = allWeight(), retVal = 'l'):

    leaderboard = [[]]
    leaderPeptide = ([], 0)
    targ = max(spectrum)
    sp = {}
    for s in spectrum:
        if sp.has_key(s):
            sp[s] += 1
        else:
            sp[s] = 1
    while len(leaderboard) > 0:
        new = []
        for p in leaderboard:
            for l in pLis:
                x = p + [l]
                new.append(x)
        order = []
        for p in new:
            m = peptideMass(p)
            if m == targ:
                s = peptideScore(sp, cyclopeptideSpectrum(p))
                # print p, s
                if s > leaderPeptide[1]:
                    leaderPeptide = (p, s)
                    lead = [p]
                elif s == leaderPeptide[1]:
                    lead.append(p)
            if m <= targ:
                order.append(p)
        leaderboard = trim(order, sp, n)
    if retVal == 'l':
        return leaderPeptide[0]
    else:
        return lead

def spectralConvolution(lis):
    a = []
    srt = sorted(lis)
    for i, x in enumerate(srt):
        for y in srt[i+1:]:
            c = y - x
            if c > 0:
                a.append(y - x)
    return a

def spectralConvolutionCount(lis):
    d = Counter([x for x in spectralConvolution(lis) if x > 56 and x < 201])
    return sorted(d.items(), key=itemgetter(1), reverse=True)

def convolutionPeptideSequencing(lis, m, n, l='l'):
    rawL = spectralConvolutionCount(lis)
    pepdides = set()
    for v in rawL[:m]:
        pepdides.add(v[0])
        val = v[1]
    for x in rawL[m:]:
        if x[1] < val:
            break
        pepdides.add(x[0])
    ans = leaderboardCyclopeptideSequencing(lis, n, pepdides, l)
    return ans


# print translation('CCUCGUACAGAAAUCAAC')
# print peptideNumber('LEADER')
# out = leaderboardCyclopeptideSequencing(lRead('cycl.txt'), 1000, retVal='s', pLis=range(57, 201))
# for lin in out:
#     print '-'.join(str(x) for x in lin),
# sp = {}
# for s in theoreticalSpectrum('VYYEVDWTMGRQIDPDEYPIAQCTRHRATILTLPDWQM'):
#     if sp.has_key(s):
#         sp[s] += 1
#     else:
#         sp[s] = 1
# print peptideScore(Counter([0, 97, 97, 129, 129, 194, 203, 226, 226, 258, 323, 323, 323, 355, 403, 452]), linearSpectrum('PEEP'))
# for out in cyclopeptideSequencing(lRead('cycl.txt')):
#     print '-'.join(str(x) for x in out),
# x = spectralConvolution([0, 137, 186, 323])
# print(x)
# out = convolutionPeptideSequencing(lRead('cycl.txt'), 20, 1100)
# # print len(out)
# # for outL in out:
# print '-'.join(str(x) for x in out),
# for x in trimT(lReadT('cycl.txt'), lRead('dna.txt'), 6):
#     print x,
# ans = sorted([0, 71, 113, 101, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416])
# tes = theoreticalSpectrum('MLAT')
# print ans
# print tes
# if ans == tes:
#     print 'yes'
# else:
#     print 'no'
# print consistent(lRead('cycl.txt'), linearSpectrum('ETC'))
print Counter(spectralConvolution(lRead("cycl.txt")))
