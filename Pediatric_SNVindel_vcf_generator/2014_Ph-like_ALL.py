#!/usr/bin/python3
usage = 'python3 '+__file__+'<2014_Ph-like_ALL_snvindel_sum> <snvindel_supptable>' 

import sys,os,json
import subprocess as sp

#function
def GET_REFALT(chro,ref,alt,pos,fa):
	p = str(int(pos) - 1)
	AddBase = list(os.popen('samtools faidx '+fa+' chr'+chro+':'+p+'-'+p))[1].strip().upper()
	return AddBase+ref,AddBase+alt,p

def SQLITEDB(sam,chr,pos,ref,alt,dna_assay,mt,tt,mn,tn,pmid,pro,vori):
	JS = {sam:{'dna_assay':dna_assay.lower(),'pmid':pmid,'project':pro,'vorigin':vori}}
	if dna_assay == 'WGS':
		TK = 'tumor_DNA_WGS_'
		NK = 'germline_DNA_WGS_'
	elif dna_assay == 'WES':
		TK = 'tumor_DNA_WES_'
		NK = 'germline_DNA_WES_'
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

PaperID = '2014_Ph-like_ALL'
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'
#Get sample name and type
SUM = open(sys.argv[1]) #sample name info
SamTy = {}
for i in SUM:
	l = i.split('\t')
	ID = l[0]
	keyId = ID.split('_')[0]
	SamTy[keyId] = [ID,l[1]]
SUM.close()

###SNVindel

#Export snvindel file
snvindelfh = open(sys.argv[2]) #snvindel table from paper
out = open(PaperID+'_SNVindel','w')
out.write('\t'.join(['Sample','CHROM','POS','Mutant_In_Tumor','Total_In_Tumor','Mutant_In_Normal','Total_In_Normal','ReferenceAllele','MutantAllele'])+'\n')

snvindelfh.readline()
for i in snvindelfh:
	l = i.split('\t')
	if l[2] in SamTy:
		SamTypeId = SamTy[l[2]][0]
		Assay = SamTy[l[2]][1]
	else:
		print('sample ID problem: '+l[2])
		continue

	chron = l[3]
	pos = l[4].split('.')[0]
	if len(l) < 16:
		MinT = l[5].split('.')[0]
		TinT = l[6].split('.')[0]
		MinN = l[7].split('.')[0]
		TinN = l[8].split('.')[0]
		REF = l[9]
		ALT= l[10]
	else:
		MinT = l[9].split('.')[0]
		TinT = l[10].split('.')[0]
		MinN = l[11].split('.')[0]
		TinN = l[12].split('.')[0]
		REF = l[13]
		ALT= l[14]
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN): #Mut allele count > total count
		MinN,TinN = TinN,MinN
	out.write('\t'.join([SamTypeId,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,Assay])+'\n')
snvindelfh.close()
out.close()

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
		sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'25207766','pcgp','somatic')
		SQLITE[ID] = sqlitejs
	else:
		if SamNam not in SQLITE[ID]:
			sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'25207766','pcgp','somatic')
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
 
