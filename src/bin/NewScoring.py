# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 00:09:12 2020

@author: reetb
"""

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


differenceList = [int(x) for x in abs(np.diff(edgelist))]
expected = [math.pow(x,2) for x in differenceList]

l = []
l.append(np.min(differenceList))
l.append(np.average(differenceList))
l.append(np.max(differenceList))
l.append(np.sum(differenceList))
l.append(np.sum(expected) / edgeCount)
l.append(np.var(differenceList))
l.append(np.std(differenceList))

print(",".join(str(i) for i in l)) 


#with open(filename + '_Gaps.txt', 'w') as f:
#    for item in differenceList:
#        f.write("%s\n" % item)

differenceList = [math.log10(x) for x in differenceList]
expected = [math.pow(x,2) for x in differenceList]

l = []
l.append(np.min(differenceList))
l.append(np.average(differenceList))
l.append(np.max(differenceList))
l.append(np.sum(differenceList))
l.append(np.sum(expected) / edgeCount)
l.append(np.var(differenceList))
l.append(np.std(differenceList))

print(",".join(str(i) for i in l))
