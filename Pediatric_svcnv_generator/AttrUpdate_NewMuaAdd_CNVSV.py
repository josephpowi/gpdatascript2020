#!/usr/bin/python3

import argparse,sys
import json

parser = argparse.ArgumentParser()
parser.add_argument('--oldFile',help="Old CNVSVLOH file")
parser.add_argument('--CNV',help="New CNVs need to be added",default=None)
parser.add_argument('--SV',help="New SVs need to be added",default=None)
parser.add_argument('--samplearray',help="sample array and pmid")
parser.add_argument('-o','--output',help="output file name",default="CNVSVLOH")
args = parser.parse_args()

if not args.oldFile or not args.samplearray:
	parser.print_help(sys.stderr)
	sys.exit(1)
#_________________________________________________________________________
#Functions
def ASSAYRETURN(attlist):
	if attlist[0]:
		return attlist[0].lower().strip()
	elif attlist[1]:
		return attlist[1].lower().strip()
	else:
		return False
def PMIDADD(JSON):
	SAMNam = JSON['sample']
	if SAMNam in SampleInfo:
		if JSON['dt'] == 10:
			dna_assay = 'wgs'
		else:
			dna_assay = ASSAYRETURN(SampleInfo[SAMNam])
		if dna_assay:
			if 'dna_assay' in JSON['mattr']:
			#if 'dna_assay' in JSON['mattr'] and JSON['mattr']['dna_assay'] == dna_assay:
				pass
			else:
				JSON['mattr']['dna_assay'] = dna_assay
		JSON['mattr']['pmid'] = SampleInfo[SAMNam][4]
		return JSON
	else:
		return False
		#print("no sample: "+SAMNam)
		#A.add(SAMNam)

def FORCNV(cline):
	clineL = cline.strip().split('\t')
	SAMNam = clineL[0]
	dna_assay = ASSAYRETURN(SampleInfo[SAMNam])
	if 'chr' in clineL[1]:
		chron = clineL[1]
	else:
		chron = 'chr'+clineL[1]
	cjs = {}
	cjs['dt'] = 4
	cjs['mattr']={'project':'pcgp'}
	cjs['mattr']['dna_assay'] = dna_assay
	cjs['mattr']['pmid'] = SampleInfo[SAMNam][4]
	cjs['mattr']['vorigin'] = 'somatic'
	cjs['sample'] = SAMNam
	cjs['value'] = float(clineL[4])
	return '\t'.join([chron,clineL[2],clineL[3],json.dumps(cjs,sort_keys=True)])+'\n'
	
def FORSV(sJS):
	SAMNam = sJS['sample']
	dna_assay = ASSAYRETURN(SampleInfo[SAMNam])
	sJS['dt'] = 5
	sJS['mattr']={'project':'pcgp'}
	sJS['mattr']['dna_assay'] = dna_assay
	sJS['mattr']['pmid'] = SampleInfo[SAMNam][4]
	sJS['mattr']['vorigin'] = 'somatic'
	if 'chrB' in sJS:
		if 'chr' not in sJS['chrB']:
			sJS['chrB'] = 'chr' + sJS['chrB']
	elif 'chrA' in sJS:
		if 'chr' not in sJS['chrA']:
			sJS['chrA'] = 'chr' + sJS['chrA']
	return json.dumps(sJS,sort_keys=True)

#pmidrnaAssay is a prebuild data structure include pmid and RNA assay type
#pmid_RNA = {}
#with open('pmidrnaAssay','rb') as f:
#	pmid_RNA = pickle.load(f)

#Extract sample information including sequencing array,project and pmid
samfh = open(args.samplearray)
SampleInfo = {}
for sam in samfh:
	samL = sam.strip().split('\t')
	if samL[0] not in SampleInfo:
		SampleInfo[samL[0]] = samL[1:]
	else:
		SampleInfo[samL[0]][4] = SampleInfo[samL[0]][4]+','+samL[5]
		for x,ap in enumerate(samL[1:5]):
			if ap:
				SampleInfo[samL[0]][x] = ap
			else:
				continue
samfh.close()

out=open(args.output,'w')
#A = set()
#Update old file
Vfh = open(args.oldFile)
for i in Vfh:
	L = i.strip().split('\t')
	js = json.loads(L[3])
	NLine = PMIDADD(js)
	if NLine:
		Newjs = json.dumps(NLine,sort_keys=True)
		out.write('\t'.join(L[0:3]+[Newjs])+'\n')
	else:
		Newjs = json.dumps(js,sort_keys=True)
		out.write('\t'.join(L[0:3]+[Newjs])+'\n')
Vfh.close()

#ADD new CNVS
CNVfh = open(args.CNV)
for c in CNVfh:
	FCNV = FORCNV(c)
	out.write(FCNV)
CNVfh.close()

#ADD new SVS
SVfh = open(args.SV)
for s in SVfh:
	sl = s.strip().split('\t')
	chron = 'chr'+sl[0]
	sjs = json.loads(sl[3])
	sjsOut = FORSV(sjs)
	out.write('\t'.join([chron,sl[1],sl[2],sjsOut])+'\n')
SVfh.close()
out.close()
