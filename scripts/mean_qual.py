import sys
import numpy as np

with open(sys.argv[1],'r') as file:
	qualities=file.read()[:-1]

a=np.array([])
for i in qualities:
	a=np.append(a, [ord(i)-33])
print(round(np.mean(a),5))
