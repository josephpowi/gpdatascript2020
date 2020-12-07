#!/usr/bin/python3
usage = 'python3 '+__file__+'<snvindel_supptable>'

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

PaperID = '2012_RB'
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'
###SNVindel

#Export snvindel file

snvindelfh = open(sys.argv[1]) #snvindel table from paper
out = open(PaperID+'_SNVindel','w')
out.write('\t'.join(['Sample','CHROM','POS','Mutant_In_Tumor','Total_In_Tumor','Mutant_In_Normal','Total_In_Normal','ReferenceAllele','MutantAllele'])+'\n')
snvindelfh.readline()

for i in snvindelfh:
	if i.startswith('Class'):
		continue
	if not i.strip():
		continue
	l = i.split('\t')
	if l[2][-2] != 'D':
		continue
	SamNam = l[2][0:-2]+'_D'
	chron = l[3]
	pos = l[4]
	MinT = l[5]
	TinT = l[6]
	MinN = l[7]
	TinN = l[8]
	REF = l[14]
	ALT = l[15]
	if TinT.isdigit() and MinT.isdigit() and int(TinT) < int(MinT): #Mut allele count > total count
		MinT,TinT = TinT,MinT
	if TinN.isdigit() and MinN.isdigit() and int(TinN) < int(MinN): #Mut allele count > total count
		MinN,TinN = TinN,MinN
	out.write('\t'.join([SamNam,chron,pos,MinT,TinT,MinN,TinN,REF,ALT,'WGS'])+'\n')
snvindelfh.close()
out.close()

#lift hg18 to hg19
fh = open(PaperID+'_SNVindel')
out = open(PaperID+'_SNVindel_hg19','w')
out.write(fh.readline().strip()+'\n')

for i in fh:
	l = i.strip().split('\t')
	tem_out = open('liftover.bed','w')
	tem_out.write(l[1]+'\t'+l[2]+'\t'+str(int(l[2])+2)+'\n')
	tem_out.close()
	sp.run("liftOver "+'liftover.bed hg18ToHg19.over.chain.gz '+'liftover_liftout liftover_nofound >tem 2>&1',shell=True)
	if os.path.getsize('liftover_nofound') == 0:
		os.remove('liftover_nofound')
		os.remove('liftover.bed')
		lifted = open('liftover_liftout')
		lifted_list = lifted.readline().strip().split('\t')
		chron = lifted_list[0]
		pos = lifted_list[1]
		out.write('\t'.join([l[0],chron,pos]+l[3:])+'\n')
		lifted.close()
		os.remove('liftover_liftout')
fh.close()
out.close()
os.rename(PaperID+'_SNVindel_hg19',PaperID+'_SNVindel')

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
		sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'22237022','pcgp','somatic')
		SQLITE[ID] = sqlitejs
	else:
		if SamNam not in SQLITE[ID]:
			sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'22237022','pcgp','somatic')
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

