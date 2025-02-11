# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:01:47 2020

@author: reetb
"""
import sys
import numpy as np
import math

filename = sys.argv[1]

edgelist = np.genfromtxt(filename, delimiter=' ')

vertexCount = int(max(max(edgelist,key=lambda item:item[1])[0], max(edgelist,key=lambda item:item[1])[1]))
edgeCount = len(edgelist)

print('*******************************************')
print('General Graph: Characteristics :')
print('*******************************************')

print('Number of vertices   :  ', vertexCount)
print('Number of edges      :  ', edgeCount)
print('*******************************************')

differenceList = [int(x) for x in abs(np.diff(edgelist))]
expected = [math.pow(x,2) for x in differenceList]

print('Minimum gap          :  ', "{0:.6f}".format(np.min(differenceList)))
print('Average gap          :  ', "{0:.6f}".format(np.average(differenceList)))
print('Maximum gap          :  ', "{0:.6f}".format(np.max(differenceList)))
print('Total (sum) gap score:  ', "{0:.6f}".format(np.sum(differenceList)))
print('Expected value of X^2:  ', "{0:.6f}".format(np.sum(expected) / edgeCount))
print('Variance is          :  ', "{0:.6f}".format(np.var(differenceList)))
print('Standard deviation   :  ', "{0:.6f}".format(np.std(differenceList)))
print('*******************************************')

with open(filename + '_Gaps.txt', 'w') as f:
    for item in differenceList:
        f.write("%s\n" % item)

differenceList = [math.log10(x) for x in differenceList]
expected = [math.pow(x,2) for x in differenceList]

print('Minimum log gap      :  ', "{0:.6f}".format(np.min(differenceList)))
print('Average log gap      :  ', "{0:.6f}".format(np.average(differenceList)))
print('Maximum log gap      :  ', "{0:.6f}".format(np.max(differenceList)))
print('Total (sum) gap score:  ', "{0:.6f}".format(np.sum(differenceList)))
print('Expected value of X^2:  ', "{0:.6f}".format(np.sum(expected) / edgeCount))
print('Variance is          :  ', "{0:.6f}".format(np.var(differenceList)))
print('Standard deviation   :  ', "{0:.6f}".format(np.std(differenceList)))
print('*******************************************')
