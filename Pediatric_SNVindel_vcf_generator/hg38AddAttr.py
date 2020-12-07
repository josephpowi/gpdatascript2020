#!/usr/bin/python3

import json,sys

#Smple table and dnaassay
SAMASSAY = {}
PUMID = {}
samfh = open(sys.argv[1])
for line in samfh:
	lineL = line.strip().split('\t')
	if lineL[0] not in SAMASSAY:
		SAMASSAY[lineL[0]] = set()
	if lineL[1].strip():
		SAMASSAY[lineL[0]].add(lineL[1].strip().lower())
	if lineL[2].strip():
		SAMASSAY[lineL[0]].add(lineL[2].strip().lower())
	if not SAMASSAY[lineL[0]]:
		SAMASSAY.pop(lineL[0])
	if lineL[0] not in PUMID:
		PUMID[lineL[0]] = lineL[5]
	else:
		if lineL[5] not in PUMID[lineL[0]]:
			PUMID[lineL[0]] += ','+ lineL[5]
samfh.close()



#add dna_assay and pmid for hg38 snvindels
vcfh = open(sys.argv[2])
out = open(sys.argv[2]+'.new','w')
for i in vcfh:
	l = i.strip().split('\t')
	ID = '@'.join(l[0:4])
	JS = json.loads(l[4])
	JSCP = JS.copy()
	for sam in JSCP:
		JS[sam]['project'] = 'pcgp'
		JS[sam]['vorigin'] = 'somatic'
		if sam not in SAMASSAY: #Only use samples that are associated with DNA assay and PMID
			JS.pop(sam)
			continue
		else:
			JS[sam]['pmid'] = PUMID[sam]
			if 'wgs' in SAMASSAY[sam]:
				assay = 'wgs'
			elif 'wes' in SAMASSAY[sam]:
				assay = 'wes'
			elif 'cc' in SAMASSAY[sam]:
				assay = 'cc'
			else:
				print('Error: not valid dna_assay info for: '+sam,file=sys.stderr)
				sys.exit(1)
			JS[sam]['dna_assay'] = assay
		if assay == 'wgs':
			TK = 'tumor_DNA_WGS_'
			NK = 'germline_DNA_WGS_'
		elif assay == 'wes':
			TK = 'tumor_DNA_WES_'
			NK = 'germline_DNA_WES_'
		elif assay == 'cc':
			TK = 'tumor_DNA_CC_'
			NK = 'germline_DNA_CC_'
		if 'tumor_DNA_ref' in JS[sam]:
			JS[sam][TK+'ref'] = JS[sam]['tumor_DNA_ref']
			JS[sam].pop('tumor_DNA_ref')
		if 'tumor_DNA_alt' in JS[sam]:
			JS[sam][TK+'alt'] = JS[sam]['tumor_DNA_alt']
			JS[sam].pop('tumor_DNA_alt')
		if 'germline_DNA_ref' in JS[sam]:
			JS[sam][NK+'ref'] = JS[sam]['germline_DNA_ref']
			JS[sam].pop('germline_DNA_ref')
		if 'germline_DNA_alt' in JS[sam]:
			JS[sam][NK+'alt'] = JS[sam]['germline_DNA_alt']
			JS[sam].pop('germline_DNA_alt')
	l[4] = json.dumps(JS)
	out.write('\t'.join(l)+'\n')
out.close()
vcfh.close()
			

