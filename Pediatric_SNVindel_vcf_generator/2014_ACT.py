#!/usr/bin/python3
usage = 'python3 '+__file__+'<snvindel_supptable3a> <snvindel_supptable3b> <snvindel_supptable3c>'

import sys,os,json

#function
def GETRA(oricha):
        orichal= oricha.split('!')[0]
        if orichal == '>':
                return ''
        else:
                return orichal

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

PaperID = '2014_ACT'
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'
###SNVindel

#Export snvindel file

#s3a
snvindelfh3a = open(sys.argv[1]) #snvindel3a table from paper
out = open(PaperID+'_SNVindel','w')
out.write('\t'.join(['Sample','CHROM','POS','Mutant_In_Tumor','Total_In_Tumor','Mutant_In_Normal','Total_In_Normal','ReferenceAllele','MutantAllele'])+'\n')

snvindelfh3a.readline()
for i in snvindelfh3a:
	l = i.split('\t')
	SamTypeId = l[0]
	chron = l[3].split('.')[0]
	pos = l[4].split('.')[0]
	MinT = ''
	TinT = ''
	MinN = ''
	TinN = ''
	REF = l[8]
	ALT= l[9]
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN): #Mut allele count > total count
		MinN,TinN = TinN,MinN
	out.write('\t'.join([SamTypeId,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,'WGS'])+'\n')
snvindelfh3a.close()
	

#s3b	
snvindelfh3b = open(sys.argv[2]) #snvindel3b table from paper
snvindelfh3b.readline()

for i in snvindelfh3b:
	l = i.split('\t')
	SamTypeId = l[0].split('.')[0]
	chron = l[3].split('.')[0]
	pos = l[4].split('.')[0]
	MinT = ''
	TinT = ''
	MinN = ''
	TinN = ''
	REF = l[8]
	ALT= l[9]
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN): #Mut allele count > total count
		MinN,TinN = TinN,MinN
	out.write('\t'.join([SamTypeId,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,'WES'])+'\n')

snvindelfh3b.close()

#s3c
snvindelfh3c = open(sys.argv[3]) #snvindel3c table from paper
snvindelfh3c.readline()
indelassay = 'WGS'

for i in snvindelfh3c:
	l = i.split('\t')
	SamTypeId = l[0].split('!')[0]
	if SamTypeId == 'SJACT009_D':
		indelassay = 'WES'
	chron = l[4].split('!')[0]
	pos = l[5].split('!')[0]
	MinT = ''
	TinT = ''
	MinN = ''
	TinN = ''
	REF = GETRA(l[9])
	ALT= GETRA(l[10])
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN): #Mut allele count > total count
		MinN,TinN = TinN,MinN
	out.write('\t'.join([SamTypeId,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,indelassay])+'\n')
snvindelfh3c.close()
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
		sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'25743702','pcgp','somatic')
		SQLITE[ID] = sqlitejs
	else:
		if SamNam not in SQLITE[ID]:
			sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'25743702','pcgp','somatic')
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

