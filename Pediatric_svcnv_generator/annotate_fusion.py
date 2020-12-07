#!/usr/bin/python3

import sys,json,pickle

sampleAssay = open(sys.argv[1])
fusionfh = open(sys.argv[2])
newfusion = open(sys.argv[3])
out = open(sys.argv[4],'w')

RNAAssay = {'23583981': 'polyA+', '22237106': 'polyA+', '23153540': 'polyA+', '25730765': 'polyA+', '24705251': 'polyA+', '24553141': 'polyA+', '25207766': 'unknown', '24332040': 'polyA+', '24703847': 'unknown', '25743702': 'polyA+', '25268584': 'unknown'}

def ASSAYRETURN(attlist):
	pmidL = []
	assayL = set()
	for pid in attlist[4].split(','):
		pmidL.append(pid)
		if pid in RNAAssay:
			assayL.add(RNAAssay[pid])
	if len(assayL) == 1:
		return list(assayL)[0]
	elif len(assayL) > 1:
		for m in assayL:
			if 'polyA' in m:
				return m
	elif len(assayL) < 1:
		return False
def FORFU(fJS):
	SAMNam = fJS['sample']
	if 'vorigin' not in fJS['mattr']:
		fJS['mattr']['vorigin'] = 'somatic' 
	if SAMNam in SampleInfo:
		rna_assay = ASSAYRETURN(SampleInfo[SAMNam])
		if rna_assay:
			fJS['mattr']['rna_assay'] = rna_assay
		if "pmid" not in fJS['mattr']:
			fJS['mattr']['pmid'] = SampleInfo[SAMNam][4]
	return json.dumps(fJS,sort_keys=True)
		
def FORNEWFU(nfJS):
	SAMNam = nfJS['sample']
	rna_assay = ASSAYRETURN(SampleInfo[SAMNam])
	del nfJS['isfusion']
	nfJS['dt'] = 2
	nfJS['mattr']={'project':'pcgp'}
	if rna_assay:
		nfJS['mattr']['rna_assay'] = rna_assay
	nfJS['mattr']['pmid'] = SampleInfo[SAMNam][4]
	nfJS['mattr']['vorigin'] = 'somatic'
	if 'chrB' in nfJS:
		if 'chr' not in nfJS['chrB']:
			nfJS['chrB'] = 'chr'+nfJS['chrB']
	elif 'chrA' in nfJS:
		if 'chr' not in nfJS['chrA']:
			nfJS['chrA'] = 'chr'+nfJS['chrA']
	return json.dumps(nfJS,sort_keys=True)

SampleInfo = {}
for sam in sampleAssay:
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
sampleAssay.close()
#update old fusion data

for f in fusionfh:
	fl = f.strip().split('\t')
	fjs = json.loads(fl[3])
	fusion = FORFU(fjs)
	out.write('\t'.join(fl[0:3]+[fusion])+'\n')

fusionfh.close()


#Add new fusion data
for nf in newfusion:
	nfl = nf.strip().split('\t')
	nfjs = json.loads(nfl[3])
	nfusion = FORNEWFU(nfjs)
	if 'chr' not in nfl[0]:
		chron = 'chr'+nfl[0]
	else:
		chron = nfl[0]
	out.write('\t'.join([chron,nfl[1],nfl[2],nfusion])+'\n')
newfusion.close()
out.close()
