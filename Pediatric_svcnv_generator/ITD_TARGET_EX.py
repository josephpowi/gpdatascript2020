#!/usr/bin/python3
"""ITD_TARGET_EX.py: usage [python3 ITD_TARGET_EX.py ITD_table TARGET_ITD(output)]
"""


import argparse,sys,re,json

parser = argparse.ArgumentParser()
parser.add_argument("-f","--file",help = "Table file containing ITD info from PanTarget table")
parser.add_argument("--extfile",help="ITD for FLT3, NOTCH1, MYC")
parser.add_argument("-m","--masterfile",help = "TARGET master file")
parser.add_argument("-o", "--output", help = "Output file name")
parser.add_argument("--DiaSample",help="pantarget diagnosis only sample list")
args = parser.parse_args()

#____________________________________________________________
#functions
def SET_POS(a,b):
	NumA = int(a.split('.')[0])
	NumB = int(b.split('.')[0])
	if NumA < NumB:
		return str(NumA),str(NumB)
	else:
		return str(NumB),str(NumA)

def GET_AALEN(aa):
	if re.search(">\d+aa",aa):
		return int(re.search(">(\d+)aa",aa).group(1))
	else:
		return len(aa.split('>')[1])

#____________________________________________________________
#Get sample info
SAMPLE = {}
SamFh = open(args.masterfile)
SamFh.readline()
for s in SamFh:
	sl = s.strip().split('\t')
	SAMPLE[sl[7]+'@'+sl[6]] = sl[5]
SamFh.close()

#Get all diagnosis samples
DIASAM = []
dsamfh = open(args.DiaSample)
for i in dsamfh:
	l = i.strip().split(' ')
	DIASAM.append(l[1])
dsamfh.close()
	




fh = open(args.file)#fusion table in pediatric.hg19.db
OUT = open(args.output,'w')
ALLITD = []
for f in fh:
	fl = f.strip().split('\t')
	chron = 'chr'+fl[5].split('.')[0]
	SamName = SAMPLE[fl[4]+'@'+fl[3]]
	if SamName not in DIASAM:
		continue
	posA, posB = SET_POS(fl[6],fl[8])
	aalen = GET_AALEN(fl[13])
	js = {}
	js['sample'] = SamName
	js['gene'] = fl[11]
	js['isoform'] = fl[14]
	js['aaduplength'] = aalen
	js['dt'] = 6
	js['mattr'] = {"project": "pantarget", "dna_assay": "cgi", "vorigin": "somatic","pmid":"29489755"}
	if '@'.join([chron,posA,posB,SamName]) in ALLITD:
		continue
	ALLITD.append('@'.join([chron,posA,posB,SamName]))
	OUT.write('\t'.join([chron,posA,posB,json.dumps(js,sort_keys=True)])+'\n')

fh.close()

#ITD for FLT3, NOTCH1, MYC
itd2fh = open(args.extfile)
for itd in itd2fh:
	itdl = itd.strip().split('\t')
	js = json.loads(itdl[3])
	if js['sample'] not in DIASAM:
		continue
	js['mattr'] = {"project": "pantarget", "dna_assay": "cgi", "vorigin": "somatic","pmid":"29489755"}
	if '@'.join(itdl[0:3]+[js['sample']]) in ALLITD:
		continue
	ALLITD.append('@'.join(itdl[0:3]+[js['sample']]))
	OUT.write('\t'.join(itdl[0:3] + [json.dumps(js,sort_keys=True)])+'\n')
itd2fh.close()
OUT.close()

