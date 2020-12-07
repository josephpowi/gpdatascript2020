#!/usr/bin/python3

import json,sys,os
import re

usage='python3 db2SQLite.py <sample_assay_table> <.tsv file dumped from sqlite3 db> <output file>'

if len(sys.argv) !=4:
	print(usage)
	sys.exit(1)
#functions
def IDENTI_REFALT(all1,all2,all1tf,all2tf,all1rnaC,all2rnaC):
	if all1tf == 't':
		ref = all1.replace('-','').strip()
		rrnaC = all1rnaC.strip()
		alt = all2.replace('-','').strip()
		mrnaC = all2rnaC.strip()
	elif all2tf == 't':
		ref = all2.replace('-','').strip()
		rrnaC = all2rnaC.strip()
		alt = all1.replace('-','').strip()
		mrnaC = all1rnaC.strip()
	else:
		print(all1,all2,all1tf,all2tf,all1rnaC,all2rnaC)
		sys.exit(1)
	return ref,alt,rrnaC,mrnaC

def GET_REFALT(chro,ref,alt,pos,fa):
	p = str(int(pos) - 1)
	AddBase = list(os.popen('samtools faidx '+fa+' chr'+chro+':'+p+'-'+p))[1].strip().upper()
	return AddBase+ref,AddBase+alt,p
def SQLITEDB(sam,chr,pos,ref,alt,dna_assay,mt,tt,mn,tn,pmid,pro,vori,rrnat,mrnat):
	JS = {sam:{'dna_assay':dna_assay.lower(),'pmid':pmid,'project':pro,'vorigin':vori}}
	if dna_assay == 'WGS':
		TK = 'tumor_DNA_WGS_'
		NK = 'germline_DNA_WGS_'
	elif dna_assay == 'WES':
		TK = 'tumor_DNA_WES_'
		NK = 'germline_DNA_WES_'
	elif dna_assay == 'CC':
		TK = 'tumor_DNA_CC_'
		NK = 'germline_DNA_CC_'
	else:
		print("dna assay not WGS no WES: "+dna_assay)
	if mt and tt:
		mt = mt.replace('-','')
		tt = tt.replace('-','')
		rt = int(tt) - int(mt)
		JS[sam][TK+'ref'] = rt
		JS[sam][TK+'alt'] = int(mt)
	if mn and tn:
		mn = mn.replace('-','')
		tn = tn.replace('-','')
		rn = int(tn) - int(mn)
		JS[sam][NK+'ref'] = rn
		JS[sam][NK+'alt'] = int(mn)
	if rrnat and mrnat:
		rrnat = rrnat.replace('-','')
		mrnat = mrnat.replace('-','')
		JS[sam]['tumor_RNA_ref'] = int(rrnat)
		JS[sam]['tumor_RNA_alt'] = int(mrnat)
	return JS


RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'

#sample not found 
#unvalidSample =['SJHYPO086_D', 'SJLGG001277_D1', 'SJLGG010124_D1', 'SJHYPO096_D', 'SJHYPO128_D', 'SJHYPO083_D', 'SJHYPO061_D', 'SJLGG001259_D1', 'SJLGG001259_D1', 'SJLGG001273_D1', 'SJHYPO112_D', 'SJHYPO065_D', 'SJHYPO065_D', 'SJHYPO093_D', 'SJHYPO080_D', 'SJHYPO101_D', 'SJLGG010119_D1', 'SJHYPO102_D', 'SJHYPO075_D', 'SJLGG001284_D1', 'SJLGG010132_D1', 'SJLGG010115_D1', 'SJHYPO103_D', 'SJHYPO049_D', 'SJHYPO091_D', 'SJHYPO122_D', 'SJHYPO023_D', 'SJHYPO064_D', 'SJLGG001271_D1', 'SJHYPO068_D', 'SJHYPO062_D', 'SJHYPO060_D', 'SJHYPO081_D', 'SJHYPO017_D', 'SJLGG010116_D1', 'SJHYPO084_D', 'SJHYPO033_D', 'SJLGG010126_D1', 'SJHYPO104_D', 'SJHYPO077_D', 'SJLGG010123_D1', 'SJLGG010117_D1', 'SJHYPO108_D', 'SJHYPO078_D', 'SJHYPO059_D', 'SJLGG001267_D1', 'SJLGG010133_D1', 'SJHYPO092_D', 'SJHYPO097_D', 'SJLGG001281_D1', 'SJLGG010131_D1', 'SJHYPO087_D', 'SJHYPO127_D', 'SJHYPO070_D', 'SJHYPO053_D', 'SJHYPO043_D', 'SJHYPO090_D', 'SJLGG001269_D1', 'SJLGG001266_D1', 'SJHYPO112_D', 'SJTALL175_D', 'SJTALL021_D', 'SJTALL163_D', 'SJTALL014_D', 'SJTALL014_D', 'SJTALL205_D', 'SJTALL205_D', 'SJTALL209_D', 'SJTALL210_D', 'SJTALL204_D', 'SJTALL162_D', 'SJTALL179_D', 'SJTALL161_D', 'SJTALL194_D', 'SJTALL195_D', 'SJTALL188_D', 'SJTALL202_D', 'SJTALL198_D', 'SJTALL199_D', 'SJBALL020383_D1', 'SJBALL021081_D1', 'SJBALL021299_D1', 'SJBALL021324_D1', 'SJBALL021333_D1', 'SJBALL021336_D1', 'SJBALL021372_D1', 'SJBALL021405_D1', 'SJBALL021427_D1', 'SJBALL021492_D1', 'SJBALL021511_D1', 'SJBALL021512_D1', 'SJBALL021514_D1', 'SJBALL087_D', 'SJBALL199_D', 'SJBALL215_D', 'SJBALL227_D', 'SJBALL230_D', 'SJBALL241_D', 'SJEWS010417_D1', 'SJEWS010418_D1', 'SJEWS010419_D1', 'SJEWS010422_D1', 'SJEWS010424_D1', 'SJEWS010425_D1', 'SJACT030_D', 'SJINF038_D', 'SJINF044_D', 'SJINF045_D', 'SJINF046_D', 'SJINF054_D', 'SJINF054_D', 'SJNBL076_D', 'SJEWS010416_D1', 'SJOS001101_M1', 'SJOS001101_M2', 'SJOS001101_M3','SJOS001101_M8','SJOS001101_M6','SJOS001101_M5','SJOS001101_M7','SJOS001101_M4','SJMB078_D','SJRHB057_R','SJEWS010423_D1','SJEWS010423_D1','SJRB001130_M1','SJTALL165_D','SJEWS010420_D1','SJRHB014_D']
unvalidSample = []

#Sample assay
SamAssay = {}
safh = open(sys.argv[1]) #sample assay table
for i in safh:
	l = i.strip().split('\t')
	if l[0] not in SamAssay:
		SamAssay[l[0]] = set()
	assay1 = l[1].strip()
	assay2 = l[2].strip()
	if assay1:
		SamAssay[l[0]].add(assay1)
	if assay2:
		SamAssay[l[0]].add(assay2)
safh.close()	

fh = open(sys.argv[2]) #.tsv file dumped from sqlite3 db
fh.readline()
SQLITE = {}
sample_noassay = set()

for line in fh:
	line = line.replace('\n','')
	l = line.split('\t')
	if l[58].upper() != 'PCGP':
		continue
	if l[7].lower() != 'somatic':
		continue
	SamNam = l[2]
	seqPlatform = l[37].strip()
	if SamNam in unvalidSample:
		continue
	if seqPlatform:
		if seqPlatform in ['{NEXT_GEN_FREQCAP}','{NEXT_GEN_EXCAP}','{NEXT_GEN_FREQEXCAP}','{NEXT_GEN_WES}']:
			ASSAY = 'WES'
		elif seqPlatform in ['{NEXT_GEN_CAPTURE}','{SANGER}']:
			ASSAY = 'CC'
		elif seqPlatform == '{NEXT_GEN_WGS}':
			ASSAY = 'WGS'
		elif l[34].strip() in ['SANGER_INDEL_DETECT','SANGER_SNP_DETECT']:
			ASSAY = 'CC'
	elif SamNam in SamAssay:
		if 'WGS' in SamAssay[SamNam]:
			ASSAY = 'WGS'
		elif 'WES' in SamAssay[SamNam]:
			ASSAY = 'WES'
		elif 'CC' in SamAssay[SamNam]:
			ASSAY = 'CC'
		#else:
		#	print('No assay for sample: '+SamNam,file=sys.stderr)
		#	sys.exit(1)
			#ASSAY = 'WGS'
	else:
		sample_noassay.add(SamNam)
		continue
		#print('sample '+SamNam+' has no dna_assay type!',file=sys.stderr)
		#sys.exit(1)
	if l[28].startswith('chr'):
		chron = l[28][3:]
	else:
		chron = l[28]
	if l[44] == l[45] and l[46] == l[47]: ###
		continue
	RefTem,AltTem,rinRNA,minRNA = IDENTI_REFALT(l[44],l[45],l[46],l[47],l[59],l[60])
	PosTem = l[29]
	if not RefTem or not AltTem:
		REF,ALT,POS = GET_REFALT(chron,RefTem,AltTem,PosTem,RefGenomeFa)
	else:
		REF,ALT,POS = RefTem,AltTem,PosTem
	ID = '@'.join([chron,POS,REF,ALT])
	if ID == '2@223163303@A@C' and SamNam == 'SJRHB046_D':
		continue
	if int(l[48]) > int(l[49]):
		#print(l[48],l[49],l[50],l[51])
		continue
	sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,ASSAY,l[48],l[49],l[50],l[51],l[57],'pcgp','somatic',rinRNA,minRNA)
	if ID not in SQLITE:
		SQLITE[ID] = sqlitejs
	else:
		if SamNam not in SQLITE[ID]:
			SQLITE[ID][SamNam] = sqlitejs[SamNam]
		else:
			newjs = sqlitejs[SamNam].copy()
			if ASSAY == 'WGS' and SQLITE[ID][SamNam]['dna_assay'] == 'WGS':
				Cur_MinT = int(l[49])
				CONT = json.dumps(SQLITE[ID][SamNam])
				print(CONT)
				MC = int(re.search('tumor_DNA_\w+?_alt":\s?(\d+)',CONT).group(1))
				RC = int(re.search('tumor_DNA_\w+?_ref":\s?(\d+)',CONT).group(1))
				if MC+RC > Cur_MinT:
					SQLITE[ID].pop(SamNam)
					sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,ASSAY,l[48],l[49],l[50],l[51],l[57],'pcgp','somatic',rinRNA,minRNA)
					SQLITE[ID][SamNam] = sqlitejs[SamNam]
				print(SQLITE[ID])
				print(SamNam)
				print(ID,ASSAY)
				print(l[48],l[49],l[50],l[51])
				#############
				###correct read count for SNVs from WGS platform and have two records from sqlite database
				#{'SJCBF004_D': {'tumor_RNA_ref': 74, 'germline_DNA_WGS_ref': 199, 'dna_assay': 'wgs', 'germline_DNA_WGS_alt': 1, 'tumor_DNA_WGS_ref': 242, 'tumor_RNA_alt': 27, 'project': 'pcgp', 'pmid': '24710217', 'vorigin': 'somatic', 'tumor_DNA_WGS_alt': 90}}
				#SJCBF004_D
				#11@85988145@T@C WGS
				#7 32 0 48 ###select small ones
				#############
			else:
				if ASSAY == 'WGS':
					SQLITE[ID][SamNam].pop('dna_assay')
				elif SQLITE[ID][SamNam]['dna_assay'] == 'WGS':
					newjs.pop('dna_assay')
				elif ASSAY == 'WES' and not SQLITE[ID][SamNam]['dna_assay'] in ['WGS','WES']:
					SQLITE[ID][SamNam].pop('dna_assay')
				else:
					newjs.pop('dna_assay') 
				if 'pmid' in newjs and 'pmid' in SQLITE[ID][SamNam] and not SQLITE[ID][SamNam]['pmid'] == newjs['pmid']:
					newpmid = ','.join([SQLITE[ID][SamNam]['pmid'],newjs['pmid']])
					SQLITE[ID][SamNam]['pmid'] = newpmid
					newjs.pop('pmid')
				SQLITE[ID][SamNam].update(newjs)
				#continue
fh.close()

out = open(sys.argv[3],'w')
for id in SQLITE:
	chron,POS,REF,ALT = id.split('@')
	out.write('\t'.join([chron,POS,REF,ALT,json.dumps(SQLITE[id],sort_keys=True)])+'\n')
out.close()

print('samples without dna_assay type:')
print(sample_noassay)
