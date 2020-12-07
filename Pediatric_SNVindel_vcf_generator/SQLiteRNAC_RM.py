#!/usr/bin/python3

import sys,json

#function
def RNACRM(jsc):
	JS = json.loads(jsc)
	JSC = {}
	for sam in JS:
		samjs = JS[sam]
		#if sam in INCLUSAMS: #
		#	JSC[sam] = samjs#
		#else:#
		if 'tumor_RNA_alt' in samjs:
			samjs.pop('tumor_RNA_alt')
		if 'tumor_RNA_ref' in samjs:
			samjs.pop('tumor_RNA_ref')
		JSC[sam] = samjs
	return JSC


##TEM include SAMPLE
#fh = open('/research/rgs01/resgen/legacy/gb_customTracks/tp/hg19/TARGET/DNA/snp.ase/rna/vcf_problemsample')
#INCLUSAMS = [x.strip() for x in fh]
#fh.close()
#FAILSAM = ['SJAML040506_R1', 'SJAML040538_R1', 'SJAML040576_R1', 'SJAML040586_R1', 'SJAML040615_R1', 'SJAML040636_R1', 'SJAML040667_R1', 'SJBALL021170_D2', 'SJCOGALL010859_R2', 'SJCOGALL010877_R2', 'SJCOGALL010904_R2', 'SJCOGALL010914_R2', 'SJNBL017366_R1']
#for s in FAILSAM:
#	if s in INCLUSAMS:
#		INCLUSAMS.remove(s)

fh = open(sys.argv[1])
out = open(sys.argv[1]+'_new','w')



for i in fh:
	i = i.replace('\n','')
	l = i.split('\t')
	if len(l[2]) != 1 or len(l[3]) != 1:
		out.write(i+'\n')
		continue
	NEWJS = RNACRM(l[4])
	l[4] = json.dumps(NEWJS)
	out.write('\t'.join(l)+'\n')
out.close()
fh.close()
