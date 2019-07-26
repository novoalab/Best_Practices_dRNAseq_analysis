#!/usr/bin/env python
import sys
import numpy as np
from scipy.stats import mannwhitneyu

def calculateErrorProbs():
    '''
    Build dictionary of error probabilities from possible Phred scores
    '''
    errorProbs = {}
    for i in range (33,74):
        errorProbs[chr(i)] = 10 ** ((33-i)/10.)
    return errorProbs
    '''
    data=line.split()
    # 0:ref
    # 1:pos
    # 2:refBase
    # 3:depth
    # 4:read(s) base(s)
    # 5:base(s) qualities
    '''
readbases = {}
qualicontainer = {}
pileup1 = sys.argv[1]

print("meanQual\tmedianQual")
with open (pileup1,'r') as p: # amplified sample 
	for line in p:
		line.rstrip("\n")
		data = line.split ()
		key = "\t".join(data[0:3])
		qualicontainer[key] = data[5].split (",")
		ary1 = np.array ([int(x) for x in qualicontainer[key]])
		mean1 = np.mean(ary1)
		median1 = np.median (ary1)
	#print key + "\t" + readbases[key]+"\t" + ",".join (qualicontainer[key]) +"\t", 
	#print data[4] + "\t" + data[5] + "\t" ,
		print("\t".join([str(x) for x in [mean1,median1]]))


        

