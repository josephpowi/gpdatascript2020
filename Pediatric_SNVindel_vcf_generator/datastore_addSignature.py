#!/usr/bin/python3

import argparse,json,sys
import subprocess as sp

parser = argparse.ArgumentParser(description='Add mutation signature to SNVindel data store')
parser.add_argument('-d','--data_store',help='data store file')
parser.add_argument('-s','--signature',help='signature file')
parser.add_argument('--panSigSample',help='panTarget samples with valid signatures detected') #PARIKF-diagnosis        1       3       4       8       10
parser.add_argument('--sigStregth',help='panTarget samples with stregth for each signature')
parser.add_argument('-g','--genome',help='hg19.gz genome file')
parser.add_argument('-o','--output',help='output file')
args = parser.parse_args()

#functions
def OPPBASEKEY(t,subt):
	BasePair = {'A':'T','T':'A','C':'G','G':'C'}
	T = t.split('>')
	newt = '>'.join([BasePair[x] for x in T])
	SUBT = list(reversed(list(subt)))
	newsubt = ''.join([BasePair[x] for x in SUBT])
	return '@'.join([newt,newsubt])
#return signature for each text in specific sample
def GETSIG(id,opid,idl): 
	samIdSig = {}
	for sam in SAMSTRE:
		newid = id+'@'+sam
		newopid = opid+'@'+sam
		sMp = [n*SAMSTRE[sam][x] for x,n in enumerate(idl)]
		sig = str(sMp.index(sorted(sMp)[-1])+1)
		samIdSig[newid] = sig
		samIdSig[newopid] = sig
	return samIdSig

def GET_BASE(chron,range):
	BASES = sp.run('samtools faidx '+RefGenome+' ' + chron+':'+range,shell=True,stdout=sp.PIPE).stdout.decode('utf-8').strip().split('\n')[1]
	return BASES.upper()

def INDEL_ADD_SIG(js,sig):
	outjs = js.copy()
	for sam in js:
		if sam in TarSigSamples:
			outjs[sam]['pantargetsignature'] = sig
	return outjs

def SNV_ADD_SIG(id,js):
	outjs = js.copy()
	for sam in js:
		if sam in TarSigSamples:
			nid = id+'@'+sam
			sig = SigNat[nid]
			if sig in TarSigSamples[sam]:
				outjs[sam]['pantargetsignature'] = sig
			else:
				outjs[sam]['pantargetsignature'] = 'n'
	return outjs


RefGenome = args.genome

#pantarget samples with valid signature 
#PARIKF-diagnosis        1       3       4       8       10
tsamfh = open(args.panSigSample)
TarSigSamples = {}
for i in tsamfh:
	l = i.strip().split('\t')
	TarSigSamples[l[0]] = l[1:]
tsamfh.close()

#panTarget samples with stregth for each signature
#CAAABF-diagnosis        0.0169003955411722      0       0       0.151384394102841       0.213951815893563       0       0       0       0       0       0
strefh = open(args.sigStregth)
SAMSTRE = {}
for s in strefh:
	l = s.strip().split('\t')
	SAMSTRE[l[0]] = list(map(float,l[1:]))
strefh.close()


#readin signature
sfh = open(args.signature)
Head = sfh.readline().strip().split('\t')
SigNat = {} #ID@SAMPLE : Signature
for line in sfh:
	L = line.strip().split('\t')
	ID = '@'.join(L[0:2])
	OPPID = OPPBASEKEY(L[0],L[1])
	PL = list(map(float,L[14:25]))
	SIGDIC = GETSIG(ID,OPPID,PL)
	SigNat.update(SIGDIC)
sfh.close()
#print(SigNat)

#Add mutation signature
fh = open(args.data_store)
out = open(args.output,'w')
for i in fh:
	I = i.replace('\n','')
	l = I.split('\t')
	JS = json.loads(l[4])
	if len(l[2].strip()) >1 or len(l[3].strip()) > 1:
		NEWINDELJS = INDEL_ADD_SIG(JS,'n')
		l[4] = json.dumps(NEWINDELJS)
		out.write('\t'.join(l)+'\n')
		continue
	chron = 'chr'+l[0]
	posRange = str(int(l[1]) - 1) + '-' + str(int(l[1]) + 1)
	mutSubtype = GET_BASE(chron,posRange)
	mutType = l[2]+'>'+l[3]
	SID = '@'.join([mutType,mutSubtype])
	NEWJS = SNV_ADD_SIG(SID,JS)
	l[4] = json.dumps(NEWJS)
	out.write('\t'.join(l)+'\n')
fh.close()
out.close()
