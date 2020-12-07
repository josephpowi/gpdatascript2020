#!/usr/bin/python3
usage = 'python3 '+__file__+' <snvindel_supptable> '+'<2011_ETP_TALL_SAMPLE>'

import sys,os,json
import subprocess as sp

#function
def GET_REFALT(chro,ref,alt,pos,fa):
	p = str(int(pos) - 1)
	AddBase = list(os.popen('samtools faidx '+fa+' chr'+chro+':'+p+'-'+p))[1].strip().upper()
	return AddBase+ref,AddBase+alt,p

#def SQLITEDB(sam,chr,pos,ref,alt,dna_assay,mt,rt,mn,rn,pmid,pro,vori):
#	JS = {sam:{'dna_assay':dna_assay,'pmid':pmid,'project':pro,'vorigin':vori,'tumor_DNA_WGS_ref':rt,'tumor_DNA_WGS_alt':mt,'germline_DNA_WGS_ref':rn,'germline_DNA_WGS_alt':mn}}
#	return '\t'.join([chr,pos,ref,alt,json.dumps(JS,sort_keys=True)])

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


#dna assay
DASSAY = {}
STFH = open(sys.argv[2])
for stl in STFH:
	STL = stl.split('\t')
	if STL[1]:
		DASSAY[STL[0]] = STL[1]
	else:
		DASSAY[STL[0]] = STL[2]
STFH.close()

PaperID = '2011_ETP_TALL'
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'
###SNVindel

#Export snvindel file

snvindelfh = open(sys.argv[1]) #snvindel table from paper
out = open(PaperID+'_SNVindel','w')
out.write('\t'.join(['Sample','CHROM','POS','Mutant_In_Tumor','Total_In_Tumor','Mutant_In_Normal','Total_In_Normal','ReferenceAllele','MutantAllele'])+'\n')
snvindelfh.readline()

for i in snvindelfh:
	l = i.split('\t')
	chron = l[3]
	pos = l[4].split('.')[0]
	MinT = l[9].split('.')[0]
	TinT = l[10].split('.')[0]
	MinN = l[11].split('.')[0]
	TinN = l[12].split('.')[0]
	REF = l[18]
	ALT= l[19]
	if int(TinT) < int(MinT):
		MinT,TinT = TinT,MinT
	if int(TinN) < int(MinN):
		MinN,TinN = TinN,MinN
	out.write('\t'.join([l[2]+'_D',chron,pos,MinT,TinT,MinN,TinN,REF,ALT,DASSAY[l[2]+'_D']])+'\n')
snvindelfh.close()
out.close()

#lift hg18 to hg19
fh = open(PaperID+'_SNVindel')
out = open(PaperID+'_SNVindel_hg19','w')
out.write(fh.readline().strip()+'\n')

for i in fh:
	l = i.strip().split('\t')
	#print(l)
	tem_out = open('liftover.bed','w')
	tem_out.write(l[1]+'\t'+l[2]+'\t'+str(int(l[2])+2)+'\n')
	tem_out.close()
	sp.run("liftOver "+'liftover.bed hg18ToHg19.over.chain.gz '+'liftover_liftout liftover_nofound >tem 2>&1',shell=True)
	if os.path.getsize('liftover_nofound') == 0:
		os.remove('liftover_nofound')
		os.remove('liftover.bed')
		lifted = open('liftover_liftout')
		lifted_list = lifted.readline().strip().split('\t')
		#print(lifted_list)
		chron = lifted_list[0]
		pos = lifted_list[1]
		out.write('\t'.join([l[0],chron,pos]+l[3:])+'\n')
		lifted.close()
		os.remove('liftover_liftout')
fh.close()
out.close()
os.rename(PaperID+'_SNVindel_hg19',PaperID+'_SNVindel')

#Generate SQLiteVCF DB
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
		sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'22237106','pcgp','somatic')
		SQLITE[ID] = sqlitejs
	else:
		if SamNam not in SQLITE[ID]:
			sqlitejs = SQLITEDB(SamNam,chron,POS,REF,ALT,l[9],l[3],l[4],l[5],l[6],'22237106','pcgp','somatic')
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
