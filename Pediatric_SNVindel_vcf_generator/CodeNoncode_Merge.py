#!/usr/bin/python3

import json,sys
import argparse

parser = argparse.ArgumentParser(description='merge pantarget coding and noncoding SNVindel data store files')
parser.add_argument('-c','--code',help='pantarget coding SNVindel data store file')
parser.add_argument('-n','--noncode',help='pantarget noncoding SNVindel data store file')
parser.add_argument('-o','--output',help='output file')
args = parser.parse_args()

#########################
#function
def ADDMATTR(js,DAssay):
	rtjs = {}
	for sam in js:
		#if 'dna_assay' in js[sam]:
		#	break
		temjs = js[sam].copy()
		if 'dna_assay' not in temjs:
			if DAssay == 'cgi':
				temjs['dna_assay'] = 'cgi'
			elif DAssay == 'wes':
				temjs['dna_assay'] = 'wes'
		temjs['project'] = 'pantarget'
		temjs['vorigin'] = 'somatic'
		temjs['pmid'] = '29489755'
		rtjs[sam] = temjs
	if not rtjs:
		#print(js)
		rtjs = js
	return rtjs

def MERGE(SQLJS,CURJS):
	sqljs = SQLJS.copy()
	curjs = CURJS.copy()
	for sam in curjs:
		if sam not in sqljs:
			sqljs[sam] = curjs[sam]
		else:
			newsqljs = COMPARE_EACHKEY(sqljs[sam],curjs[sam])
			sqljs[sam] = newsqljs
	return sqljs

def COMPARE_EACHKEY(S1JS,S2JS):
	s1js = S1JS.copy()
	s2js = S2JS.copy()
	for k in s2js:
		if k not in s1js:
			s1js[k] = s2js[k]
	return s1js


##data merge
SQLITE = {}
#Noncoding
NFH = open(args.noncode)
for i in NFH:
	l = i.strip().split('\t')
	ID = '@'.join(l[0:4])
	JS = json.loads(l[4])
	NJS = ADDMATTR(JS,'cgi')
	SQLITE[ID] = NJS
NFH.close()

#Coding
CFH = open(args.code)
for line in CFH:
	L = line.strip().split('\t')
	ID = '@'.join(L[0:4])
	JS = json.loads(L[4])
	CJS = ADDMATTR(JS,None)
	if ID not in SQLITE:
		SQLITE[ID] = CJS
	else:
		MERJS = MERGE(SQLITE[ID],CJS)
		SQLITE[ID] = MERJS
CFH.close()

#output merge file
out = open(args.output,'w')
for id in SQLITE:
	if not SQLITE[id]:
		continue
	IDL = id.split('@')
	OUTL = IDL+[json.dumps(SQLITE[id])]
	out.write('\t'.join(OUTL)+'\n')
out.close()
