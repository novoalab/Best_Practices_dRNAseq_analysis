#!/usr/bin/env python
import sys
# ~ import numpy as np
# ~ from scipy.stats import mannwhitneyu

# ~ def mann_whitney_test (list1,list2):
    # ~ ary1 = np.array(list1)
    # ~ ary2 = np.array(list2)
    # ~ return mannwhitneyu(ary1,ary2)

# ~ def cal_man_whitney_z_score (samp1, samp2):
    # ~ s1_len = samp1.shape[0]
    # ~ s2_len = samp2.shape[0]
    # ~ tot_len = s1_len + s2_len

    # ~ all_vals = np.concatenate([samp1,samp2])
    # ~ ranks = np.empty (tot_len, int)
    # ~ ranks[all_vals.argsort()] = np.arange(1,tot_len+1)
    # ~ s1_ranks_sum = ranks[:s1_len].sum()
    # ~ u1 = s1_ranks_sum - (s1_len * (s1_len + 1)) / 2

    # ~ mu = s1_len * s2_len / 2
    # ~ rhou =  np.sqrt (s1_len * s2_len * (s1_len + s2_len + 1)/12)
# ~ #    z = np.abs(u1 - mu) /

# ~ bases = {'A','G','C','T'}
# ~ pileupbases = bases.union({',','.'})
# ~ dec_quals = {}

# ~ def calculateErrorProbs():
    # ~ '''
    # ~ Build dictionary of error probabilities from possible Phred scores
    # ~ '''
    # ~ errorProbs = {}
    # ~ for i in range (33,74):
        # ~ errorProbs[chr(i)] = 10 ** ((33-i)/10.)
    # ~ return errorProbs


def process(line):
    '''
    calls snp(s) for a single line
    '''
    data=line.split()
    # 0:ref
    # 1:pos
    # 2:refBase
    # 3:depth
    # 4:read(s) base(s)
    # 5:base(s) qualities

    #refbase = data[2].upper()
    qualitites = data[5]
    newline = '\t'.join(data[0:5])
    # iterate through bases from mapped reads (a pileup column)
    qualiIdx = 0
    dec_q = ""
    for q in data[5]:
        dec_q += str (ord(q)-33) + ','
    dec_q = dec_q.rstrip (",")
    print(newline+"\t"+dec_q)	#I have added the parenthesis


pileup = sys.argv[1]
with open (pileup,'r') as p:
	for line in p:
		line = line.rstrip("\n")
		data = line.split()
		process(line)			#I have added a tab here
#        if  int(data[3]) > 19:
#            process (line)

