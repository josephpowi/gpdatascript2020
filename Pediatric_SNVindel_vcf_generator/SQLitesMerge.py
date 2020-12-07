#!/usr/bin/python3

import json,sys
import argparse

parser = argparse.ArgumentParser(description='merge SNVindel data store files')
parser.add_argument('--SQLites',nargs='+',help='SNVindel data store files need to be merged')
parser.add_argument('-o','--output',help='output merged file')
parser.add_argument('--anno',action='store_true',help='If anno is specified, vep annotation info will be placed at the 6th column!')
parser.add_argument('--rsid',action='store_true',help='If rsid is specified, SNP ID will be placed at the 7th column! rsid will automatically add "." as vepanno if --anno if not specified.')
args = parser.parse_args()

#functions
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
	if 'vorigin' in s1js and 'vorigin' in s2js and s1js['vorigin'] != s2js['vorigin']:
		print('vorigin: '+s1js['vorigin'] + ' ' +s2js['vorigin'])
		#print(S1JS,S2JS)
		sys.exit(1)
	if 'pmid' in s2js:
		if 'pmid' not in s1js:
			s1js['pmid'] = s2js['pmid']
		if s1js['pmid'] != s2js['pmid']:
			A = s1js['pmid'].split(',')
			B = s2js['pmid'].split(',')
			C = list(set(A).union(set(B)))
			F = ','.join(C)
			s1js['pmid'] = F
	if 'dna_assay' in s1js and 'dna_assay' in s2js and s1js['dna_assay'] != s2js['dna_assay']:
		if 'wgs' in s1js['dna_assay'] or 'wgs' in s2js['dna_assay']:
			s1js['dna_assay'] = 'wgs'
		elif 'cgi' in s1js['dna_assay'] or 'cgi' in s2js['dna_assay']:
			s1js['dna_assay'] = 'cgi'
		elif 'wes' in s1js['dna_assay'] or 'wes' in s2js['dna_assay']:
			s1js['dna_assay'] = 'wes'
		elif 'cc' in s1js['dna_assay'] or 'cc' in s2js['dna_assay']:
			s1js['dna_assay'] = 'cc'
		else:
			print('dna_assay error: '+s1js['dna_assay']+','+s2js['dna_assay'],file=sys.stderr)
			sys.exit(1)
	for k in s2js:
		if k not in s1js:
			s1js[k] = s2js[k]
	return s1js
	
#data merge
SQLITE = {}
if args.anno:
	VEPANNO = {}
else:
	VEPANNO = False

if args.rsid:
	RSID = {}
else:
	RSID = False

for f in args.SQLites:
	fh = open(f)
	for i in fh:
		#print(i)
		i = i.replace('\n','')
		l = i.strip().split('\t')
		ID = '@'.join(l[0:4])
		JS = json.loads(l[4])
		if isinstance(VEPANNO,dict):
			if ID not in VEPANNO or VEPANNO[ID] == '.':
				try:
					VEPANNO[ID] = l[5]
				except:
					VEPANNO[ID] = '.'
		if isinstance(RSID,dict):
			if ID not in RSID or RSID[ID] == '.':
				try:
					RSID[ID] = l[6]
				except:
					RSID[ID] = '.'
		if ID not in SQLITE:
			SQLITE[ID] = JS
			continue
		else:
			NEWJS = MERGE(SQLITE[ID],JS)
			SQLITE[ID] = NEWJS
	fh.close()


#output merged file
out = open(args.output,'w')
for id in SQLITE:
	if not SQLITE[id]:
		continue
	IDL = id.split('@')
	OUTL = IDL+[json.dumps(SQLITE[id])]
	if args.anno:
		if id not in VEPANNO:
			print(id +'not in VEP')
			sys.exit(1)
		OUTL.append(VEPANNO[id])
	if args.rsid:
		if len(OUTL) == 5:
			OUTL.append('.')
		if id in RSID:
			OUTL.append(RSID[id])
		else:
			OUTL.append('.')
	out.write('\t'.join(OUTL)+'\n')

out.close()
