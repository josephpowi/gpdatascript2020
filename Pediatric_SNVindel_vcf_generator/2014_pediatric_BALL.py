#!/usr/bin/python3
usage = 'python3 '+__file__+' <sample_match> <snvindel_supptable>'

import sys,os,json

#function
def GETINFO(barco,lineL,Tumor):
	chron = lineL[3]
	pos = lineL[4].split('.')[0]
	MinT = lineL[11].split('.')[0]
	TinT = lineL[12].split('.')[0]
	MinR = lineL[13].split('.')[0]
	TinR = lineL[14].split('.')[0]
	MinN = lineL[15].split('.')[0]
	TinN = lineL[16].split('.')[0]
	Ref = lineL[9]
	Alt = lineL[10]
	if Tumor == 'Both':
		return '\t'.join([sam_match[barco][0],chron,pos,MinT,TinT,MinN,TinN,Ref,Alt,'WES']) + '\n' + '\t'.join([sam_match[barco][1],chron,pos,MinR,TinR,MinN,TinN,Ref,Alt,'WES'])
	elif Tumor == 'D':
		return '\t'.join([sam_match[barco][0],chron,pos,MinT,TinT,MinN,TinN,Ref,Alt,'WES'])
	elif Tumor == 'R':
		return '\t'.join([sam_match[barco][1],chron,pos,MinR,TinR,MinN,TinN,Ref,Alt,'WES'])

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

PaperID = '2014_pediatric_BALL'
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'
###SNVindel

sam_match = {}
samMatfh = open(sys.argv[1])
for i in samMatfh:
	l = i.strip().split('\t')
	if l[0] not in sam_match:
		sam_match[l[0]] = ['D','R']
	if '_D1' in l[1]:
		sam_match[l[0]][0] = l[1]
	elif '_R1' in l[1]:
		sam_match[l[0]][1] = l[1]
samMatfh.close()


#Export snvindel file
snvindelfh = open(sys.argv[2]) #snvindel table from paper
out = open(PaperID+'_SNVindel','w')
out.write('\t'.join(['Sample','CHROM','POS','Mutant_In_Tumor','Total_In_Tumor','Mutant_In_Normal','Total_In_Normal','ReferenceAllele','MutantAllele'])+'\n')
snvindelfh.readline()


for i in snvindelfh:
	l = i.split('\t')
	OUTINFO = GETINFO(l[2],l,l[0])
	out.write(OUTINFO)
	out.write('\n')
out.close()
snvindelfh.close()


#output sqlite db file
fh = open(PaperID+'_SNVindel')
SQLITE = {}
for i in fh:
	if i.startswith('Sample'):
		continue
	i = i.replace('\n','')
	l = i.split('\t')
	SamNam = l[0]
	if SamNam.endswith('_D1'):
		continue
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
		sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'25790293','pcgp','somatic')
		SQLITE[ID] = sqlitejs
	else:
		if SamNam not in SQLITE[ID]:
			sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'25790293','pcgp','somatic')
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
