#!/usr/bin/python3

"""not useful...
TARGET_NONCOD_MAF2VCF.py: Usage-(python3 TARGET_NONCOD_MAF2VCF.py TARGET_SNVindel_NONCODE.vcf VCFHEAD targetid2st ~/tp/anno/db/target.samples /research/rgs01/project_space/zhanggrp/PanTARGET/common/yanling/pan_target/CGI/annovar_rescue_filterStep3/deliver.all_tier23 /research/rgs01/project_space/zhanggrp/PanTARGET/common/yanling/pan_target/CGI/annovar_rescue_filterStep3/deliver.aml_tier23 /research/rgs01/project_space/zhanggrp/PanTARGET/common/yanling/pan_target/CGI/annovar_rescue_filterStep3/deliver.nbl_tier23 /research/rgs01/project_space/zhanggrp/PanTARGET/common/yanling/pan_target/CGI/annovar_rescue_filterStep3/deliver.os_tier23 /research/rgs01/project_space/zhanggrp/PanTARGET/common/yanling/pan_target/CGI/annovar_rescue_filterStep3/deliver.wt_tier23)
                             This script takes as many maf files and combines and trasfers to vcf file.
                             The first args is the output file.
                             The second args is VCF head file.
                             The third args is targetid2st.
                             The Fourth args is target.sample file  followed by all your maf file.
"""


import sys,re,os,json,argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c','--CGIID',help='CGI ID from Diagnosis samples')
parser.add_argument('--sampleinfo',help='sample infomations with 6 digits id and SJid')
parser.add_argument('-o','--output',help='output file name')
parser.add_argument('files',nargs='+',help='all snvindel files from tier23')
args = parser.parse_args()



####################################
#function
def GET_REFALT(chro,ref,alt,pos,fa):
	if not ref or not alt:
		p = str(int(pos) - 1)
		AddBase = list(os.popen('samtools faidx '+fa+' chr'+chro+':'+p+'-'+p))[1].strip().upper()
		return AddBase+ref,AddBase+alt,p
	else:
		return ref,alt,pos


#Decide which one is alternative allele
def GET_RefAlt(ReAl,TuAl1,TuAl2):
	if TuAl1 == ReAl:
		return ReAl.replace('-','').strip(),TuAl2.replace('-','').strip()
	else:
		return ReAl.replace('-','').strip(),TuAl1.replace('-','').strip()



#Extraction of read count for each svnindel in both normal and tumor samples (json format)
def GET_DEPTH(tf,ta,nr,na):
	DEP = {}
	DEP['tumor_DNA_CGI_ref'] = int(tf)
	DEP['tumor_DNA_CGI_alt'] = int(ta)
	DEP['germline_DNA_CGI_ref'] = int(nr)
	DEP['germline_DNA_CGI_alt'] = int(na)
	return DEP

#reference genome sequence hg19
RefGenomeFa = '/research/rgs01/resgen/legacy/gb_customTracks/tp/genomes/hg19.gz'

####################################
#all diagnosis CGIID
cidfh = open(args.CGIID)
DIAGCGID = []
for cid in cidfh:
	DIAGCGID.append(cid.strip())
cidfh.close()


#Sample master file (target.samples) info
#6 digits ID vs SJID
TID2SJID = {}
SAMINFO = open(args.sampleinfo)
for i in SAMINFO:
	l = i.strip().split(' ')
	TID2SJID[l[0]] = l[1]
SAMINFO.close()

####################################
#readin maf files
#Generate data store file
InFh = args.files
SNVindelDic = {}
for maf in InFh:
	MAF = open(maf)
	for m in MAF:
		ml = m.split('\t')
		if len(ml) < 51:
			m = m.rstrip() + MAF.readline()
			ml = m.split('\t')
		TSamBar = ml[31]
		if TSamBar not in DIAGCGID:
			continue
		D6ID = TSamBar.split('-')[2]
		if D6ID in TID2SJID:
			SamId =  TID2SJID[D6ID]
		else:
			#print(D6ID)
			continue
		temREF,temALT = GET_RefAlt(ml[26].strip(),ml[27].strip(),ml[28].strip())
		chron,tempos = ml[34],ml[35]
		REF,ALT,pos = GET_REFALT(chron,temREF,temALT,tempos,RefGenomeFa)
		if pos in ['61975631','61975630','61975632']:
			print('1'+REF+'1','1'+ALT+'1',pos,chron)
		TuCountRef,TuCountAlt,NorCountRef,NorCountAlt = ml[19],ml[18],ml[22],ml[21]
		CurVar = '\t'.join([chron,pos,REF,ALT])
		SamJS = GET_DEPTH(TuCountRef,TuCountAlt,NorCountRef,NorCountAlt)
		if SamJS:
			if CurVar not in SNVindelDic:
				SNVindelDic[CurVar] = {SamId:SamJS}
			else:
				SNVindelDic[CurVar][SamId] = SamJS
		else:
			continue
	MAF.close()


###################################################
#OUTPUT
OUT = open(args.output,'w')
for var in SNVindelDic:
	OUT.write('\t'.join([var,json.dumps(SNVindelDic[var])])+'\n')
OUT.close()

