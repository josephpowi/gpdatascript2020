#!/usr/bin/python3
"""parse_cnv.py: USAGE: python3 collect_cnv.py original_file<SJBALL020016_C1_Single_CONSERTING_NoCREST_Mapability_100.txt>
		   output file is hard coded as CNV
"""

import sys,json,os


if len(sys.argv) <2:
	os.system("sed -n '2,4'p parse_cnv.py")
	sys.exit(1)

#________________________________________
#Function
def GET_COR(p1,p2):
	if int(p1) < int(p2):
		return p1,p2
	elif int(p1) > int(p2):
		return p2,p1

oriCNV = open(sys.argv[1])
OUT = open("cnv",'w')

oriCNV.readline()
for line in oriCNV:
	lineL = line.strip().split('\t')
	chron = 'chr'+lineL[0]
	posS,posE = GET_COR(lineL[1],lineL[2])
	log2ratio = lineL[5]
	if log2ratio == '0':
		continue
	else:
		js = {}
		js['dt'] = 4
		js['sample'] = 'SJHYPO010000_C213-Nalm16'
		js['value'] = float(log2ratio)
		js['mattr'] = {"project":"pedccl"}
		js['mattr']['dna_assay'] = 'wgs'
		js['mattr']['vorigin'] = 'somatic'
		jsout = json.dumps(js,sort_keys=True)
		OUT.write('\t'.join([chron,posS,posE,jsout])+'\n')
oriCNV.close()
OUT.close()
