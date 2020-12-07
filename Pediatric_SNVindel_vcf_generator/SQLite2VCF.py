#!/usr/bin/python3

import argparse,json,sys

parser = argparse.ArgumentParser(description="Convert SNVindel data store to vcf file")
parser.add_argument('-j','--JFfile',help='json format SNVindel file (SQLite)')
parser.add_argument('--header',help='prebuilt vcf header file')
parser.add_argument('-o','--output',help='output vcf file')
parser.add_argument('--anno',action='store_true',help='If anno is specified, vep annotation info will be extracted from the 6th column!')
parser.add_argument('--rsid',action='store_true',help='If rsid is specified, SNP ID will be extracted from the 7th column!')
args=parser.parse_args()

#functions
def PARSEJS(samvaljs):
	SamVal = []
	GETSAMPLE(samvaljs.keys()) #function GETSAMPLE
	PARSEJS_FORMAT_List = SETFORMAT(samvaljs)
	PARSEJS_FORMAT = ':'.join(PARSEJS_FORMAT_List)
	SamVal.append(PARSEJS_FORMAT)
	for s in samvaljs:
		samvallist = []
		for attrkey in PARSEJS_FORMAT_List:
			if 'tumor' in attrkey or 'germline' in attrkey:
				if attrkey+'_ref' in samvaljs[s] and attrkey+'_alt' in samvaljs[s]:
					samvallist.append(str(samvaljs[s][attrkey+'_ref'])+','+str(samvaljs[s][attrkey+'_alt']))
				elif attrkey+'_ref' in samvaljs[s] or attrkey+'_alt' in samvaljs[s]:
					print('only have reads count for ref or alt: ',file=sys.stderr)
				else:
					samvallist.append('.')
			#elif attrkey == 'mutation_signature':
			#	if attrkey in samvaljs[s]:
			#		sigjs = samvaljs[s][attrkey]
			#		samk = list(sigjs.keys())[0]
			#		SIGFormatV = sigjs[samk]
			#		samvallist.append('|'.join([samk,SIGFormatV]))
			#	else:
			#		samvallist.append('.')
			else:
				if attrkey in samvaljs[s]:
					samvallist.append(samvaljs[s][attrkey])
				else:
					samvallist.append('.')
		if len(SamVal) == 1:
			SamVal.append([s])
			SamVal.append([':'.join(samvallist)])
		elif len(SamVal) > 1:
			SamVal[1].append(s)
			SamVal[2].append(':'.join(samvallist))
		else:
			print('no format generated for:',file=sys.stderr)
	#print(SamVal)
	return SamVal
		

#get final format for vcf file 
def SETFORMAT(fjs):
	setformat_format = ['tumor_DNA_WGS',\
			'germline_DNA_WGS',\
			'tumor_DNA_CGI',\
			'germline_DNA_CGI',\
			'tumor_DNA_WES',\
			'germline_DNA_WES',\
			'tumor_DNA',\
			'germline_DNA',\
			'tumor_DNA_CC',\
			'germline_DNA_CC',\
			'tumor_RNA',\
			'dna_assay',\
			'pantargetsignature',\
			'tcgaskcmsignature',\
			'project',\
			'vorigin',\
			'pmid']
	allkey = set()
	for sam in fjs.keys():
		for k in fjs[sam]:
			if 'tumor' in k or 'germline' in k:
				testk = k[:-4]
			else:
				testk = k
			allkey.add(testk)
	format_list = []
	for fk in setformat_format:
		if fk in allkey:
			format_list.append(fk)
	return format_list
		

def GETSAMPLE(sams):
	for s in sams:
		if s in SAMPLES:
			continue
		else:
			SAMPLES.append(s)


#output file
out = open(args.output,'w')


#header file
headerfh = open(args.header)
for i in headerfh:
	out.write(i.strip()+'\n')
headerfh.close()

#if vep annotation provided from SQLite file
#VEPANNO="""##VEP="v88" time="2018-06-10 14:33:16" cache="/home/jwang7/.vep/homo_sapiens_refseq/91_GRCh37" ensembl-variation=88.1e40b4b ensembl-funcgen=88.905998a ensembl-io=88.277fe7c ensembl=88.b8ff470 1000genomes="phase3" COSMIC="81" ClinVar="201706" ESP="20141103" HGMD-PUBLIC="20164" assembly="GRCh37.p13" dbSNP="150" gencode="GENCODE 19" genebuild="2011-04" gnomAD="170228" polyphen="2.2.2" refseq="01_2015" regbuild="1.0" sift="sift5.2.2"
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|HGVS_OFFSET">
#"""
VEPANNO="""##VEP="v99" time="2020-08-17 20:02:02" cache="/research/rgs01/resgen/legacy/gb_customTracks/tp/jwang/tools/VEP/.vep/homo_sapiens_refseq/96_GRCh37" ensembl-io=99.441b05b ensembl-variation=99.a7f8736 ensembl-funcgen=99.0832337 ensembl=99.d3e7d31 1000genomes="phase3" COSMIC="86" ClinVar="201810" ESP="20141103" HGMD-PUBLIC="20174" assembly="GRCh37.p13" dbSNP="151" gencode="GENCODE 19" genebuild="2011-04" gnomAD="r2.1" polyphen="2.2.2" refseq="01_2015" regbuild="1.0" sift="sift5.2.2"
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|GIVEN_REF|USED_REF|BAM_EDIT|HGVS_OFFSET">
"""
if args.anno:
	out.write(VEPANNO)

#main process
SAMPLES = []
VARSamVal = {}
if args.anno:
	VEPAN = {}
if args.rsid:
	SNPDB = {}
fh = open(args.JFfile)
for line in fh:
	line = line.replace('\n','')
	lineL = line.split('\t')
	JS = json.loads(lineL[4])
	ID = '\t'.join(lineL[0:4])
	FormatSamVal = PARSEJS(JS)
	VARSamVal[ID] = FormatSamVal
	if args.anno:
		VEPAN[ID] = lineL[5]
	if args.anno and args.rsid:
		SNPDB[ID] = lineL[6]
	elif args.rsid:
		SNPDB[ID] = lineL[5]
fh.close()
print(len(SAMPLES))
#output vcf file
out.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+SAMPLES)+'\n')
for id in VARSamVal:
	if len(VARSamVal[id]) < 3:
		print(id,VARSamVal[id])
		continue
	outformat = VARSamVal[id][0]
	outsams = VARSamVal[id][1]
	outvals = VARSamVal[id][2]
	chron,pos,REF,ALT = id.split('\t')
	formatLength = outformat.split(':')
	outValList = [':'.join(['.']*len(formatLength))] * len(SAMPLES)
	for x,s in enumerate(outsams):
		IDX = SAMPLES.index(s)
		outValList[IDX] = outvals[x]
	outlist = [chron,pos,'.',REF,ALT,'.','.','.',outformat] + outValList
	if args.anno:
		outlist[7] = VEPAN[id]
	if args.rsid:
		outlist[2] = SNPDB[id]
	out.write('\t'.join(outlist)+'\n')
out.close()

