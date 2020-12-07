#!/usr/bin/python3

"""PanTarget_VCF_Pre.py: usage:"python3 PanTarget_VCF_Pre.py [target.samples(The sixth column SJBALL020336_D2)] snvindel_target"
                        Use pan_target_sample_info.update_Nov212016_OSbadAdded_update02082017_0A4I0S_back to replace target.samples.
                                python3 PanTarget_VCF_Pre.py [pan_target_sample_info.update_Nov212016_OSbadAdded_update02082017_0A4I0S_back snvindel_target] 
                        This script is used to prepare snvindel vcf file for all Target samples to be annotated by VEP. 
                         NOTE: ITD variations are excluded for this moment until we figure out how to annotate ITD. Because
                         VEP is not handling ITD well.
                         For some of variations without both CGI and WXS(e.g. nomut     nomut   TALL    diagnosis       PARASZ  1
                        47704968                        indel   TAL1    noncoding       SE_INS  NM_003189       miseq   miseq   miseq
                        miseq   miseq   miseq   miseq   miseq   newadded                PARASZ  1       47704968                        TALL
                        diagnosis       indel   TAL1    noncoding       SE_INS  NM_003189       nobam   nobam   nobam   nobam   nobam   nobam
                        nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   
                        nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam   nobam
                        nobam   nobam   nobam   nobam   miseq   miseq   miseq   miseq   nobam   -1      TARGET-10-PARASZ-09A-01D        SJALL015582_D1
                        -1      -1      TARGET-10-PARASZ-09A-01R        -1      -1      -1      nobam   nobam   nobam   nobam   nobam   nobam),
                        they are pending with ITD.
"""

import os
import json
import argparse


parser = argparse.ArgumentParser(description='Convert PanTarget coding table from Xiaotu to data store format file')
parser.add_argument('-s','--sampleinfo',help='pan-target Diagnosis sample information') #CAAABC CAAABC-diagnosis or PAMXHJ SJBALL020339_D2-PAMXHJ
parser.add_argument('-i','--snvindel',help='pan-target coding snvindels') #with two genome coordinates modified
parser.add_argument('-o','--output',help='output file name')
parser.add_argument('--TAL1',help='TAL1 promoter insertion',default=False) #PARASZ	1	47704967	A	AAC
args = parser.parse_args()


####################################
#function
def GET_REFALT(chro,ref,alt,pos,s_i_t,fa):
	if s_i_t == 'indel' and not ref or not alt:
		p = str(int(pos) - 1)
		AddBase = list(os.popen('samtools faidx '+fa+' chr'+chro+':'+p+'-'+p))[1].strip().upper()
		return AddBase+ref,AddBase+alt,p
	else:
		return ref,alt,pos

#Extraction of read count for each svnindel in both normal and tumor samples (json format)
#Rules:tumor_DNA_CGI:germline_DNA_CGI:tumor_DNA_WES:germline_DNA_WES:tumor_RNA
def GET_DEPTH(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10):
	DEP = {}
	if n1.isdigit() and n2.isdigit():
		DEP['tumor_DNA_CGI_ref'] = int(n2) - int(n1)
		DEP['tumor_DNA_CGI_alt'] = int(n1)
		DEP['germline_DNA_CGI_ref'] = int(n6) - int(n5)
		DEP['germline_DNA_CGI_alt'] = int(n5)
		DEP['dna_assay'] = 'cgi'
	if n3.isdigit() and n4.isdigit():
		DEP['tumor_DNA_WES_ref'] = int(n4) - int(n3)
		DEP['tumor_DNA_WES_alt'] = int(n3)
		DEP['germline_DNA_WES_ref'] = int(n8) - int(n7)
		DEP['germline_DNA_WES_alt'] = int(n7)
		if not 'dna_assay' in DEP:
			DEP['dna_assay'] = 'wes'
	if n9.isdigit() and n10.isdigit():
		DEP['tumor_RNA_ref'] = int(n10) - int(n9)
		DEP['tumor_RNA_alt'] = int(n9)
	return DEP

#reference genome sequence hg19
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'

####################################
#Sample master file (target.samples) info
SamMastF = open(args.sampleinfo)
SamMa = {}
for sm in SamMastF:
	sml = sm.strip().split(' ')
	SamMa[sml[0]] = sml[1]
SamMastF.close()

####################################
#Generate data store file
SNVIndelF = open(args.snvindel)
SNVIndelF.readline()
SNVindelDic = {}
for svi in SNVIndelF:
	svil = svi.split('\t')
	if svil[8].strip() == 'ITD':
		continue
	P6dig = svil[4].strip()
	if P6dig not in SamMa or svil[3].lower() != 'diagnosis':
		print('This SNVindel is not from Diagnosis sample: '+svi)
		continue
	SamNam = SamMa[P6dig]
	chron = svil[5].split('.')[0].strip()
	VarType = svil[9].strip()
	RefTem = svil[7].strip().replace('-','')
	AltTem = svil[8].strip().replace('-','')
	PosTem = svil[6].split('.')[0].strip()
	Info = GET_REFALT(chron,RefTem,AltTem,PosTem,VarType,RefGenomeFa)
	REF = Info[0]
	ALT = Info[1]
	POS = Info[2]
	#if svil[41].strip().startswith('dbsnp'):
	#	ID = svil[41].strip().split(':')[1]
	#else:
	#	ID = '.'
	MCGIT = svil[14].strip().split('.')[0]
	TCGIT = svil[15].strip().split('.')[0]
	MCGIN = svil[16].strip().split('.')[0]
	TCGIN = svil[17].strip().split('.')[0]
	MWXST = svil[18].strip().split('.')[0]
	TWXST = svil[19].strip().split('.')[0]
	MWXSN = svil[20].strip().split('.')[0]
	TWXSN = svil[21].strip().split('.')[0]
	MRNAT = svil[0].strip().split('.')[0]
	TRNAT = svil[1].strip().split('.')[0]
	SamJS = GET_DEPTH(MCGIT,TCGIT,MWXST,TWXST,MCGIN,TCGIN,MWXSN,TWXSN,MRNAT,TRNAT)
	if SamJS:
		CurVar = '\t'.join([chron,POS,REF,ALT])
		if CurVar not in SNVindelDic:
			SNVindelDic[CurVar] = {SamNam:SamJS}
		else:
			SNVindelDic[CurVar][SamNam] = SamJS
	else:
		continue
SNVIndelF.close()

##################################################
#TAL1
if args.TAL1:
	TFH = open(args.TAL1)
	TJS = {'dna_assay':'amplicon sequencing','project':'pantarget','vorigin':'somatic','pmid':'28671688'}
	for tline in TFH:
		TL = tline.strip().split('\t')
		SAMNAME = SamMa[TL[0]]
		ID = '\t'.join(TL[1:])
		if ID in SNVindelDic:
			SNVindelDic[ID][SAMNAME] = TJS
		else:
			SNVindelDic[ID] = {SAMNAME:TJS}
	TFH.close()	



###################################################
#OUTPUT
OUT = open(args.output,'w')
for var in SNVindelDic:
	OUT.write('\t'.join([var,json.dumps(SNVindelDic[var])])+'\n')
OUT.close()

