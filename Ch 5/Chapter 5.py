__author__ = 'maurice'
from math import factorial

def paths(m, n):
    return factorial(m+n)/(factorial(m)*factorial(n))

def coins(denom, n):
    m = n + 1
    ans = {i:0 for i in range(m)}
    vals = [i for i in denom if i < m]
    z = 1
    while vals != []:
        for i in vals:
            ans[i] = z
        new = set(j+c for j in vals for c in denom)
        vals = [p for p in new if p < m and ans[p] == 0]
        z += 1
    return ans


print coins([9,5,3,1],  19163)