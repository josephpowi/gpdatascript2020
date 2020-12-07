#!/usr/bin/python3
usage = 'python3 '+__file__+'<sample> <snvindel_supptable>'

import sys,os,json
import subprocess as sp

#function
def GET_REFALT(chro,ref,alt,pos,fa):
	p = str(int(pos) - 1)
	AddBase = list(os.popen('samtools faidx '+fa+' chr'+chro+':'+p+'-'+p))[1].strip().upper()
	return AddBase+ref,AddBase+alt,p

def SQLITEDB(sam,chr,pos,ref,alt,dna_assay,mt,rt,mn,rn,pmid,pro,vori):
	JS = {sam:{'dna_assay':dna_assay,'pmid':pmid,'project':pro,'vorigin':vori,'tumor_DNA_WGS_ref':rt,'tumor_DNA_WGS_alt':mt,'germline_DNA_WGS_ref':rn,'germline_DNA_WGS_alt':mn}}
	return '\t'.join([chr,pos,ref,alt,json.dumps(JS,sort_keys=True)])

def SQLITEDB(sam,chr,pos,ref,alt,dna_assay,mt,tt,mn,tn,pmid,pro,vori):
	mt = mt.replace('NA','')
	tt = tt.replace('NA','')
	mn = mn.replace('NA','')
	tn = tn.replace('NA','')
	JS = {sam:{'dna_assay':dna_assay.lower(),'pmid':pmid,'project':pro,'vorigin':vori}}
	if dna_assay.upper() == 'WGS':
		TK = 'tumor_DNA_WGS_'
		NK = 'germline_DNA_WGS_'
	elif dna_assay.upper() == 'WES':
		TK = 'tumor_DNA_WES_'
		NK = 'germline_DNA_WES_'
	else:
		print(dna_assay)
		sys.exit(1)
	if mt and tt:
		rt = int(tt) - int(mt)
		JS[sam][TK+'ref'] = rt
		JS[sam][TK+'alt'] = int(mt)
	if mn and tn:
		rn = int(tn) - int(mn)
		JS[sam][NK+'ref'] = rn
		JS[sam][NK+'alt'] = int(mn)
	return JS
	#return '\t'.join([chr,pos,ref,alt,json.dumps(JS,sort_keys=True)])



if len(sys.argv) == 1:
        print(usage)
        sys.exit(1)

PaperID = '2012_hypodiploid_ALL'
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'

#sample and assay
samfh = open(sys.argv[1]) #sample
samAssay = {}
for i in samfh:
	l = i.split('\t')
	if 'WES' in i:
		samAssay[l[0]] = 'wes'
	elif 'WGS' in i:
		samAssay[l[0]] = 'wgs'
samfh.close()

###SNVindel

#Export snvindel file
snvindelfh = open(sys.argv[2]) #snvindel table from paper
out = open(PaperID+'_SNVindel','w')
germout = open(PaperID+'_germ_SNVindel','w')

out.write('\t'.join(['Sample','CHROM','POS','Mutant_In_Tumor','Total_In_Tumor','Mutant_In_Normal','Total_In_Normal','ReferenceAllele','MutantAllele'])+'\n')
germout.write('\t'.join(['Sample','CHROM','POS','Mutant_In_Tumor','Total_In_Tumor','Mutant_In_Normal','Total_In_Normal','ReferenceAllele','MutantAllele'])+'\n')

snvindelfh.readline()

for i in snvindelfh:
	l = i.split('\t')
	samL = l[2].split('-')
	indiv = samL[0].strip()
	TYPEL = samL[1].split('and')
	ALLSAM = []
	for m in TYPEL:
		ALLSAM.append(indiv+'_'+m.strip())
	for SamTypeId in ALLSAM:
		chron = l[3]
		pos = l[4].split('.')[0]
		MinN = l[13].split('.')[0]
		TinN = l[14].split('.')[0]
		if '_D' in SamTypeId:
			MinT = l[9].split('.')[0]
			TinT = l[10].split('.')[0]
		elif '_R' in SamTypeId:
			MinT = l[11].split('.')[0]
			TinT = l[12].split('.')[0]
		REF = l[15]
		ALT= l[16]
		if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT):
			MinT,TinT = TinT,MinT
		if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN):
			MinN,TinN = TinN,MinN
		if l[18] == 'germline':
			germout.write('\t'.join([SamTypeId,chron,pos,MinT,TinT,MinN,TinN,REF,ALT])+'\n')
		elif chron == 'chr17' and pos == '7579313' and REF == 'G' and ALT == 'C' and SamTypeId == 'SJHYPO055_D':
			germout.write('\t'.join([SamTypeId,chron,pos,MinT,TinT,MinN,TinN,REF,ALT])+'\n')
		else:
			out.write('\t'.join([SamTypeId,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,samAssay[SamTypeId]])+'\n')
snvindelfh.close()
out.close()
germout.close()

#output sqlite db
fh = open(PaperID+'_SNVindel')
SQLITE = {}
for i in fh:
	if i.startswith('Sample'):
		continue
	i = i.replace('\n','')
	l = i.split('\t')
	SamNam = l[0]
	if l[1].startswith('chr'):
		chron = l[1][3:]
	else:
		chron = l[1]
	RefTem = l[7].replace('-','').strip()
	AltTem = l[8].replace('-','').strip()
	PosTem = l[2]
	if not RefTem or not AltTem:
		REF,ALT,POS = GET_REFALT(chron,RefTem,AltTem,PosTem,RefGenomeFa)
	else:
		REF,ALT,POS = RefTem,AltTem,PosTem
	ID = '@'.join([chron,POS,REF,ALT])
	if ID not in SQLITE:
		sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'23334668','pcgp','somatic')
		SQLITE[ID] = sqlitejs
	else:
		if SamNam not in SQLITE[ID]:
			sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'23334668','pcgp','somatic')
			SQLITE[ID][SamNam] = sqlitejs[SamNam]
		else:
			#print(SQLITE[ID])
			#print(SamNam)
			#print(ID,l[9])
			#print(l[3],l[4],l[5],l[6])
			continue
fh.close()

out = open(PaperID+'_SNVindel.sqlite','w')
for id in SQLITE:
	chron,POS,REF,ALT = id.split('@')
	if REF == ALT: #remove record having same reference and alternative alleles
		continue
	out.write('\t'.join([chron,POS,REF,ALT,json.dumps(SQLITE[id],sort_keys=True)])+'\n')
out.close()

