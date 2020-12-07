#!/usr/bin/python3
usage = 'python3 '+__file__+'<snv_supptable> <indel_supptable> <vali_snv_supptable> <vali_indel_supptable>'

import sys,os,json
import subprocess as sp

if len(sys.argv) == 1:
	print(usage)
	sys.exit(1)

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


PaperID = '2013_E_RHB'
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'
###SNVindel

#output SNV
snvfh = open(sys.argv[1]) #snv table from paper
out = open(PaperID+'_SNVindel','w')
out.write('\t'.join(['Sample','CHROM','POS','Mutant_In_Tumor','Total_In_Tumor','Mutant_In_Normal','Total_In_Normal','ReferenceAllele','MutantAllele'])+'\n')

snvfh.readline()
for i in snvfh:
	l = i.split('\t')
	SamNam = l[0]
	chron = l[3].split('.')[0]
	pos = l[4].split('.')[0]
	MinT = l[9].split('.')[0]
	RinT = l[10].split('.')[0]
	TinT = str(int(MinT) + int(RinT))
	MinN = ''
	TinN = ''
	REF,ALT = l[8].split('/')
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN): #Mut allele count > total count
		MinN,TinN = TinN,MinN
	out.write('\t'.join([SamNam,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,'WGS'])+'\n')
snvfh.close()

indelfh = open(sys.argv[2]) #indel table from paper

indelfh.readline()
for i in indelfh:
	l = i.split('\t')
	SamNam = l[0]
	chron = l[4].split('.')[0]
	pos = l[5].split('.')[0]
	GType = l[9].split('/')
	REF = GType[0].replace('-','')
	ALT = GType[1].replace('-','')
	MinT = l[10].split('.')[0].replace('-','')
	RinT = l[11].split('.')[0].replace('-','')
	TinT = str(int(MinT) + int(RinT))
	MinN = ''
	TinN = ''
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN): #Mut allele count > total count
		MinN,TinN = TinN,MinN
	out.write('\t'.join([SamNam,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,'WGS'])+'\n')
indelfh.close()


#validation cohort
valisnvfh = open(sys.argv[3])
valisnvfh.readline()
for i in valisnvfh:
	l = i.split('\t')
	if l[2][-1] in ['M','A','R']:
		SamNam = l[2][0:-1] +'_'+l[2][-1]
	else:
		SamNam = l[2]+'_D'
	chron = l[3]
	pos = l[4]
	MinT = l[9]
	TinT = l[10]
	MinN = l[11]
	TinN = l[12]
	REF = l[13]
	ALT = l[14]
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN): #Mut allele count > total count
		MinN,TinN = TinN,MinN
	out.write('\t'.join([SamNam,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,'WES'])+'\n')

valindelfh = open(sys.argv[4])
valindelfh.readline()

for i in valindelfh:
	l = i.split('\t')
	if l[2][-1] in ['M','A','R']:
		SamNam = l[2][0:-1] +'_'+l[2][-1]
	else: 
		SamNam = l[2]+'_D'
	chron = l[3]
	pos = l[4]
	MinT = l[9]
	TinT = l[10]
	MinN = l[11]
	TinN = l[12]
	REF = l[13]
	ALT = l[14]
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN): #Mut allele count > total count
		MinN,TinN = TinN,MinN
	out.write('\t'.join([SamNam,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,'WES'])+'\n')

out.close()



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
		sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'24332040','pcgp','somatic')
		SQLITE[ID] = sqlitejs
	else:
		if SamNam not in SQLITE[ID]:
			sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'24332040','pcgp','somatic')
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

