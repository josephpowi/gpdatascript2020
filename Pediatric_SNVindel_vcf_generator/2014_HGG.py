#!/usr/bin/python3
usage = 'python3 '+__file__+' <snvindel_SuppTable>'

import sys,os,json
import subprocess as sp

if len(sys.argv) == 1:
        print(usage)
        sys.exit(1)

#functions
def GET_REFALT(chro,ref,alt,pos,fa):
	p = str(int(pos) - 1)
	AddBase = list(os.popen('samtools faidx '+fa+' chr'+chro+':'+p+'-'+p))[1].strip().upper()
	return AddBase+ref,AddBase+alt,p

def SQLITEDB(sam,chr,pos,ref,alt,dna_assay,mt,tt,mn,tn,pmid,pro,vori):
	JS = {sam:{'dna_assay':dna_assay.lower(),'pmid':pmid,'project':pro,'vorigin':vori}}
	if dna_assay.upper() == 'WGS':
		TK = 'tumor_DNA_WGS_'
		NK = 'germline_DNA_WGS_'
	elif dna_assay.upper() == 'WES':
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


PaperID = '2014_HGG'
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'
###SNVindel

#Export snvindel file

snvindelfh = open(sys.argv[1]) #snvindel table from paper
out = open(PaperID+'_SNVindel','w')
out.write('\t'.join(['Sample','CHROM','POS','Mutant_In_Tumor','Total_In_Tumor','Mutant_In_Normal','Total_In_Normal','ReferenceAllele','MutantAllele'])+'\n')
snvindelfh.readline()

for i in snvindelfh:
	l = i.split('\t')
	SamTypeId = l[1]
	chron = l[2].split('.')[0]
	pos = l[3].split('.')[0]
	MinT = l[8].split('.')[0]
	TinT = l[9].split('.')[0]
	REF = l[10]
	ALT= l[11]
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	out.write('\t'.join([SamTypeId,chron,pos,MinT,TinT,'','',REF,ALT,'WES'])+'\n')
snvindelfh.close()

#WGS
wgssnvindelfh = open(sys.argv[2]) #WGS snvindel
wgssnvindelfh.readline()
for i in wgssnvindelfh:
	l = i.split('\t')
	SamTypeId = l[1]
	chron = l[2].split('.')[0]
	pos = l[3].split('.')[0]
	MinT = l[4].split('.')[0]
	TinT = l[5].split('.')[0]
	REF = l[6]
	ALT = l[7]
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	out.write('\t'.join([SamTypeId,chron,pos,MinT,TinT,'','',REF,ALT,'WGS'])+'\n')

wgssnvindelfh.close()

out.close()


#fix ref and alt alleles swap error
fh = open(PaperID+'_SNVindel')
out = open(PaperID+'_SNVindel.new','w')
out.write(fh.readline())
for i in fh:
	l = i.strip().split('\t')
	if len(l[7]) != 1 or len(l[8])!= 1 or l[7] == '-' or l[8] == '-':
		out.write(i)
		continue
	REF = l[7]
	ALT = l[8]
	RELREF = sp.run('samtools faidx '+RefGenomeFa+' chr'+l[1]+':'+l[2]+'-'+l[2],shell=True,stdout=sp.PIPE).stdout.decode('utf-8').strip().split('\n')[1]
	if REF.upper() != RELREF.upper():
		l[7] = ALT
		l[8] = REF
		out.write('\t'.join(l)+'\n')
	else:
		out.write(i)
fh.close()
out.close()
os.system('rm -f '+PaperID+'_SNVindel')
os.system('mv '+PaperID+'_SNVindel.new '+PaperID+'_SNVindel')



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
		sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'24705251','pcgp','somatic')
		SQLITE[ID] = sqlitejs
	else:
		if SamNam not in SQLITE[ID]:
			sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'24705251','pcgp','somatic')
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
