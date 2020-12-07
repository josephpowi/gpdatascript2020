#!/usr/bin/python3

import argparse
import json

parser = argparse.ArgumentParser(description='Add Germline DNA count')
parser.add_argument('-g','--GERMCount',help='./DNARNACOUNT/GERMCOUNT')
parser.add_argument('-s','--sqlitedatastore',help='PCGP snvindel data store')
args = parser.parse_args()

#function
def GERMADD(id,jsc):
	JS = json.loads(jsc)
	germjs = {}
	for sam in JS:
		temjs = JS[sam].copy()
		if 'germline_DNA_WGS_ref' in temjs or 'germline_DNA_WES_ref' in temjs or 'germline_DNA_CC_ref' in temjs:
			germjs[sam] = temjs
			continue
		if 'tumor_DNA_WGS_ref' in temjs and id in GERMCount['WGS'] and sam in GERMCount['WGS'][id]:
			temjs['germline_DNA_WGS_ref'] = GERMCount['WGS'][id][sam][0]
			temjs['germline_DNA_WGS_alt'] = GERMCount['WGS'][id][sam][1]
		if 'tumor_DNA_WES_ref' in temjs and id in GERMCount['WES'] and sam in GERMCount['WES'][id]:
			temjs['germline_DNA_WES_ref'] = GERMCount['WES'][id][sam][0]
			temjs['germline_DNA_WES_alt'] = GERMCount['WES'][id][sam][1]
		#print(id,sam,temjs,JS[sam])
		germjs[sam] = temjs
	return germjs

with open(args.GERMCount) as input:
	GERMCount = json.load(input)

#Add Germline DNA Count
out = open(args.sqlitedatastore+'.GERMDNACOUNT','w')
fh = open(args.sqlitedatastore)

for i in fh:
	i = i.replace('\n','')
	l = i.split('\t')
	ID = '@'.join(l[0:4])
	NEWJS = GERMADD(ID,l[4])
	l[4] = json.dumps(NEWJS)
	out.write('\t'.join(l)+'\n')
out.close()
fh.close()
