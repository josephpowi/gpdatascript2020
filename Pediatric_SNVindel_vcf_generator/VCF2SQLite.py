#!/usr/bin/python3

import argparse,json,gzip

parser = argparse.ArgumentParser(description="Convert vcf file to SNVindel data store")
parser.add_argument('-v','--vcf',help="vcf file")
parser.add_argument('-o','--output',help="output file")
parser.add_argument('--anno',action='store_true',help='If anno is specified, vep annotation info will be placed at the 6th column!')
parser.add_argument('--rsid',action='store_true',help='If rsid is specified, SNP ID will be placed at the 7th column!')
args = parser.parse_args()

#___________________________________________________
#function
def CRESQLITE(formt,ValList):
	JS = {}
	ForList = formt.split(':')
	for v,s in ValList:
		valList = v.split(':')
		FVPair = list(zip(ForList,valList))
		JS[s] = PARSEFORMAT(FVPair)
	return JS

def PARSEFORMAT(fv):
	sjs = {}
	for f,v in fv:
		if v == '.':
			continue
		if 'tumor' in f or 'germline' in f:
			ValL = list(map(int,v.split(',')))
			sjs[f+'_ref'] = ValL[0]
			sjs[f+'_alt'] = ValL[1]
		else:
			sjs[f] = v
	return sjs



#Process vcf file
VCF = args.vcf
if VCF.endswith('.gz'):
	Vfh = gzip.open(VCF)
else:
	Vfh = open(VCF)


SAMPLES = [] #ALL samples
out = open(args.output,'w')

for line in Vfh:
	if isinstance(line,bytes):
		line = line.decode('utf-8')
	if line.startswith('#CHROM'):
		SAMPLES = line.strip().split('\t')[9:]
		continue
	elif line.startswith('#'):
		continue
	L = line.strip().split('\t')
	chron = L[0]
	POS = L[1]
	REF = L[3]
	ALT = L[4]
	FORMAT = L[8]
	EMPTVAL = ':'.join(['.']*len(FORMAT.split(':')))
	VALIVALUE = [(v,SAMPLES[x]) for x,v in enumerate(L[9:]) if v != EMPTVAL]
	sqlJS = CRESQLITE(FORMAT,VALIVALUE)
	outline = '\t'.join([chron,POS,REF,ALT,json.dumps(sqlJS,sort_keys=True)])
	if args.anno:
		outline += '\t' + L[7]
	if args.rsid:
		outline += '\t' + L[2]
	out.write(outline+'\n')

out.close()
Vfh.close()
