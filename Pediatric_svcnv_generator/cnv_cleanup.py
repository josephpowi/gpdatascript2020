#!/usr/bin/python3

import argparse,sys

parser = argparse.ArgumentParser()
parser.add_argument('-f','--CNV',help="The raw CNV file (cnv_analysis_l2r_pre_2017-09-25_14_28_12.txt)")
parser.add_argument('-o','--OUTPUT',help="Output file",default="cnv.raw")
args=parser.parse_args()

if len(sys.argv) <3:
	parser.print_help()
	sys.exit(1)


#_____________________________________________________
#Functions
def EX_POS(P1,P2):
	if int(P1) > int(P2):
		return P2,P1
	elif int(P1) < int(P2):
		return P1,P2

#______________________________________________________
#varaints
Invalid_pos = 0
DUP = 0
CNV_LINE = []
out = open(args.OUTPUT,'w')
out_dup = open('CNVDUP','w')
out_invpos = open('CNVinvalpos','w')

#______________________________________________________
#raw CNV file processing
CNVfh = open(args.CNV)
SAM_ExcluChr8CNVS = ['SJBALL020877_D1','SJBALL020013_D1','SJBALL020625_D1']
out.write(CNVfh.readline())
for line in CNVfh:
	lineList = line.split('\t')
	posSTem = lineList[11]
	posETem = lineList[12]
	if abs(int(posSTem) - int(posETem)) == 1:
		Invalid_pos += 1
		out_invpos.write(line)
		continue
	else:
		PosS,PosE = EX_POS(posSTem,posETem)
		SamName = lineList[0]
		chron = lineList[10]
		val = lineList[13]
		if SamName in SAM_ExcluChr8CNVS and chron == '8': #remove CNVs from chr8 of three samples including SJBALL020877_D1, SJBALL020013_D1, SJBALL020625_D1
			continue
		CkID = '@'.join([chron,PosS,PosE,SamName,val])
		if CkID in CNV_LINE:
			DUP += 1
			out_dup.write(line)
			continue
		else:
			out.write(line)
			CNV_LINE.append(CkID)
CNVfh.close()
out.close()
out_dup.close()
out_invpos.close()
#Report the number of cnvs with invalid pos or duplicated cnvs
print(str(Invalid_pos) + " CNVs with invalid position were removed!", file=sys.stderr)
print(str(DUP) + " duplicated CNVs were removed!", file = sys.stderr)
