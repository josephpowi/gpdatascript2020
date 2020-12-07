#!/usr/bin/python3

import json,sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--SQLpaper',help='SQLite file from PCGP paper')
parser.add_argument('--SQLhg38',help='SQLite file from hg38 vcf',default=None)
parser.add_argument('--SQLSQL',help='SQLite file from pp sqlite db')
#parser.add_argument('--samtab',help='pre-prepared sample table')
parser.add_argument('-o','--output',help='output file')
args=parser.parse_args()

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
	if s1js['vorigin'] != s2js['vorigin']:
		print('vorigin: '+s1js['vorigin'] + ' ' +s2js['vorigin'])
		sys.exit(1)
	if 'pmid' in s2js:
		if s1js['pmid'] != s2js['pmid']:
			A = s1js['pmid'].split(',')
			B = s2js['pmid'].split(',')
			C = list(set(A).union(set(B)))
			F = ','.join(C)
			s1js['pmid'] = F
	if s1js['dna_assay'] != s2js['dna_assay']:
		if 'wgs' in s1js['dna_assay'] or 'wgs' in s2js['dna_assay']:
			s1js['dna_assay'] = 'wgs'
		elif 'wes' in s1js['dna_assay'] or 'wes' in s2js['dna_assay']:
			s1js['dna_assay'] = 'wes'
		elif 'cc' in s1js['dna_assay'] or 'cc' in s2js['dna_assay']:
			s1js['dna_assay'] = 'cc'
		else:
			print(s1js['dna_assay'],s2js['dna_assay'])
			sys.exit(1)
	for k in s2js:
		if k not in s1js:
			s1js[k] = s2js[k]
	return s1js
	
#Smple table and dnaassay
#SAMASSAY = {}
#PUMID = {}
#samfh = open(args.samtab)
#for line in samfh:
#	lineL = line.strip().split('\t')
#	if lineL[0] not in SAMASSAY:
#		SAMASSAY[lineL[0]] = set()
#	if lineL[1].strip():
#		SAMASSAY[lineL[0]].add(lineL[1].strip().lower())
#	if lineL[2].strip():
#		SAMASSAY[lineL[0]].add(lineL[2].strip().lower())
#	if not SAMASSAY[lineL[0]]:
#		SAMASSAY.pop(lineL[0])
#	if lineL[0] not in PUMID:
#		PUMID[lineL[0]] = lineL[5]
#	else:
#		if lineL[5] not in PUMID[lineL[0]]:
#			PUMID[lineL[0]] += ','+ lineL[5]
	
#samfh.close()

	
SQLITE = {}
NOREAD = []
#SNVindels from PCGP paper
paperfh = open(args.SQLpaper)
for i in paperfh:
	l = i.strip().split('\t')
	ID = '@'.join(l[0:4])
	JS = json.loads(l[4])
	SQLITE[ID] = JS
	for sam in JS:
		if len(JS[sam]) <= 4:
			NOREAD.append('@'.join(l[0:4]+[sam]))
paperfh.close()

#SNVindels from protein paint database
dbfh = open(args.SQLSQL)
for i in dbfh:
	l = i.strip().split('\t')
	ID = '@'.join(l[0:4]) 
	JS = json.loads(l[4])
	if ID not in SQLITE:
		SQLITE[ID] = JS
		continue
	else:
		NEWJS = MERGE(SQLITE[ID],JS)
		SQLITE[ID] = NEWJS
dbfh.close()
	
#SNVindels from hg38 vcf
if args.SQLhg38:
	vcfh = open(args.SQLhg38)
	for i in vcfh:
		l = i.strip().split('\t')
		ID = '@'.join(l[0:4])
		JS = json.loads(l[4])
		#JSCP = JS.copy()
		#for sam in JSCP:
		#	JS[sam]['project'] = 'pcgp'
		#	JS[sam]['vorigin'] = 'somatic'
		#	if sam not in SAMASSAY:
		#		JS.pop(sam)
		#		continue
				#assay = 'wgs'
				#JS[sam]['dna_assay'] = 'wgs' #if sample has not been found from PCGP paper, assign wgs for the moment!!!
				#print(sam)
		#	else:
		#		JS[sam]['pmid'] = PUMID[sam]
		#		if 'wes' in SAMASSAY[sam]:
		#			assay = 'wes'
		#		elif 'wgs' in SAMASSAY[sam]:
		#			assay = 'wgs'
		#		elif 'cc' in SAMASSAY[sam]:
		#			assay = 'cc'
		#		else:
		#			print(SAMASSAY[sam],sam)
		#			sys.exit(1)
		#		JS[sam]['dna_assay'] = assay
		#	if assay == 'wgs':
		#		TK = 'tumor_DNA_WGS_'
		#		NK = 'germline_DNA_WGS_'
		#	elif assay == 'wes':
		#		TK = 'tumor_DNA_WES_'
		#		NK = 'germline_DNA_WES_'
		#	elif assay == 'cc':
		#		TK = 'tumor_DNA_CC_'
		#		NK = 'germline_DNA_CC_'
		#	if 'tumor_DNA_ref' in JS[sam]:
		#		JS[sam][TK+'ref'] = JS[sam]['tumor_DNA_ref']
		#		JS[sam].pop('tumor_DNA_ref')
		#	if 'tumor_DNA_alt' in JS[sam]:
		#		JS[sam][TK+'alt'] = JS[sam]['tumor_DNA_alt']
		#		JS[sam].pop('tumor_DNA_alt')
		#	if 'germline_DNA_ref' in JS[sam]:
		#		JS[sam][NK+'ref'] = JS[sam]['germline_DNA_ref']
		#		JS[sam].pop('germline_DNA_ref')
		#	if 'germline_DNA_alt' in JS[sam]:
		#		JS[sam][NK+'alt'] = JS[sam]['germline_DNA_alt']
		#		JS[sam].pop('germline_DNA_alt')
		if ID not in SQLITE:
			SQLITE[ID] = JS
		else:
			NEWJS = MERGE(SQLITE[ID],JS)
			SQLITE[ID] = NEWJS

out = open(args.output,'w')
for id in SQLITE:
	if not SQLITE[id]:
		continue	
	IDL = id.split('@')
	out.write('\t'.join(IDL+[json.dumps(SQLITE[id])])+'\n')
