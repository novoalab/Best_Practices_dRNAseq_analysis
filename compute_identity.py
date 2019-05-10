import sys
from Bio import AlignIO

al = sys.argv[1] 
align = AlignIO.read(al, "fasta")

for n in range(0,len(align[0])):
	n = 0
	m = 0
	mis = 0
	d = 0
	i = 0
	j = 0
	ins = 0
	dels = 0
	previousi = False
	previousd = False
	while n<len(align[0]):
		column = align[:,n]
		if column[0] == column[1]:
			if column[0] != "-":
				previousi = False
				previousd = False
				m += 1
			else:
				if previousi == False:
					previousi = True
					ins += 1
				if previousd == False:
					previousd = True
					dels += 1
				j += 1 #gap-gap
				i += 1
				d += 1
		else:
			if column[0] == "-":
				previousd = False
				if previousi == False:
					previousi = True
					ins += 1
				i += 1
			elif column[1] == "-":
				previousi = False
				if previousd == False:
					previousd = True
					dels += 1
				d += 1
			else:
				previousi = False
				previousd = False
				mis += 1
		if n == 50:
			i50 = i
			ins50 = ins
			d50 = d
			dels50 = dels
		n += 1
length = len(align[0])-j

id = round(100*(m/length),2) # matches / (matches + mismatches + insertions + deletions)

if 'i50' in locals():
	print(length, m, mis, i, d, id, i50, d50, ins, dels, ins50, dels50, sep="\t")
else:
	print(length, m, mis, i, d, id, i, d, ins, dels, ins, dels, sep="\t")







