#!/usr/bin/python3

import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('-r','--RNACount',help='tp/hg19/TARGET/DNA/snp.ase/rna/TARGET_MDS_RNAC')
parser.add_argument('-s','--sqlitedatastore',help='TARGET snvindel data store')
args = parser.parse_args()

#function
def RNACADD(id,jsc):
	if id in TARGETRNAC:
		EXARC = True
	else:
		EXARC = False
	JS = json.loads(jsc)
	rnajs = {} 
	for sam in JS:
		temjs = JS[sam].copy()
		if EXARC and sam in TARGETRNAC[id]:
			temjs['tumor_RNA_ref'] = TARGETRNAC[id][sam][0]
			temjs['tumor_RNA_alt'] = TARGETRNAC[id][sam][1]
		else:
			if sam in RNASAMPLE:
				temjs['tumor_RNA_ref'] = 0
				temjs['tumor_RNA_alt'] = 0
				#print(id,sam)
				#print(temjs)
		rnajs[sam] = temjs
	return rnajs
			

with open(args.RNACount) as input:
	TARGETRNAC = json.load(input)
RNASAMPLE = []
for id in TARGETRNAC:
	CURSAM = [x for x in list(TARGETRNAC[id].keys()) if x not in RNASAMPLE]
	RNASAMPLE.extend(CURSAM)
print(len(RNASAMPLE))

#Add rna count
out = open(args.sqlitedatastore+'.RNACOUNT','w')
fh = open(args.sqlitedatastore)

for i in fh:
	i = i.replace('\n','')
	l = i.split('\t')
	if len(l[2]) != 1 or len(l[3]) != 1:
		out.write(i+'\n')
		continue
	ID = '@'.join(l[0:4])
	NEWJS = RNACADD(ID,l[4])
	l[4] = json.dumps(NEWJS)
	out.write('\t'.join(l)+'\n')
out.close()
fh.close()

