#!/usr/bin/python3

import json,sys

file = open(sys.argv[1]) #Fusion file
out = open(sys.argv[2],'w') #output file


#________________________________________________________
#Function
def GETFUGENAB(gA,gB,fujson):
	fujson['geneA'] = gA
	fujson['geneB'] = gB
	return fujson

def HANDLE_DUX4IGH(sam,pid,jsFile):
	fjsA = {}
	fjsA['sample'] = sam
	fjsA['dt'] = 2
	fjsA['mattr'] = {'project':'pcgp'}
	fjsA['mattr']['pmid'] = pid
	fjsA['mattr']['vorigin'] = 'somatic'
	fjsA['strandA'] = jsFile[0]['a']['strand']
	fjsA['strandB'] = jsFile[0]['b']['strand']
	fjsB = fjsA.copy()
	posA = (jsFile[0]['a']['codon'] -1)*3 + 190998973
	posB = 106659897
	fjsA['posA'] = posA
	fjsB['posB'] = posB
	fjsB['chrB'] = jsFile[0]['b']['chr']
	fjsA['chrA'] = 'chr4'
	if 'name' in jsFile[0]['a']:
		fugeneA = jsFile[0]['a']['name']
	else:
		fugeneA = fjsA['chrA']
	if 'name' in jsFile[0]['b']:
		fugeneB = jsFile[0]['b']['name']
	else:
		fugeneB = fjsB['chrB'] = jsFile[0]['b']['chr']
	fjsA = GETFUGENAB(fugeneA,fugeneB,fjsA)
	fjsB = GETFUGENAB(fugeneA,fugeneB,fjsB)
	#fjsA['fusiongene'] = fugeneA + ' > ' + fugeneB
	#fjsB['fusiongene'] = fugeneA + ' > ' + fugeneB
	return '\t'.join([fjsB['chrB'],str(posB),str(posB),json.dumps(fjsA,sort_keys=True)])+'\n' + '\t'.join([fjsA['chrA'],str(posA),str(posA),json.dumps(fjsB,sort_keys=True)])+'\n'

def HANDLE_TY9(sam,pid,jsFile):
	fjsA = {}
	fjsA['dt'] = 2
	fjsA['sample'] = sam
	fjsA['mattr'] = {'project':'pcgp'}
	if pid:
		fjsA['mattr']['pmid'] = str(pid)
	fjsA['mattr']['vorigin'] = 'somatic'
	fjsA['strandA'] = jsFile['strand']
	fjsA['strandB'] = jsFile['partner']['strand']
	if 'chimericreads' in jsFile:
		fjsA['clipreadA'] = jsFile['chimericreads'] if jsFile['chimericreads'] > 0 else 1
		if jsFile['ratio'] == 0:
			RatioA = 0.01
		else:
			RatioA = jsFile['ratio']
		fjsA['totalreadA'] = int((fjsA['clipreadA']) // RatioA)
	if 'chimericreads' in jsFile['partner']:
		fjsA['clipreadB'] = jsFile['partner']['chimericreads'] if jsFile['partner']['chimericreads'] > 0 else 1
		if jsFile['partner']['ratio'] == 0:
			RatioB = 0.01
		else:
			RatioB = jsFile['partner']['ratio']
		fjsA['totalreadB'] = int((fjsA['clipreadB']) // RatioB)
	fjsB = fjsA.copy()
	posA = jsFile['position']
	posB = jsFile['partner']['position']
	fjsA['posA'] = posA
	fjsB['posB'] = posB
	fjsA['chrA'] = jsFile['chr']
	fjsB['chrB'] = jsFile['partner']['chr']
	if 'gene' in jsFile:
		fugeneA = jsFile['gene']
	else:
		fugeneA = fjsA['chrA']
	if 'gene' in jsFile['partner']:
		fugeneB = jsFile['partner']['gene']
	else:
		fugeneB = fjsB['chrB']
	fjsA = GETFUGENAB(fugeneA,fugeneB,fjsA)
	fjsB = GETFUGENAB(fugeneA,fugeneB,fjsB)
	#fjsA['fusiongene'] = fugeneA + ' > ' + fugeneB
	#fjsB['fusiongene'] = fugeneA + ' > ' + fugeneB
	return '\t'.join([fjsB['chrB'],str(posB),str(posB),json.dumps(fjsA,sort_keys=True)]) + '\n' + '\t'.join([fjsA['chrA'],str(posA),str(posA),json.dumps(fjsB,sort_keys=True)])+'\n'
def HANDLE_TY8(sam,pid,jsFile):
	fjsA = {}
	fjsA['dt'] = 2
	fjsA['sample'] = sam
	fjsA['mattr'] = {'project':'pcgp'}
	if pid:
		fjsA['mattr']['pmid'] = str(pid)
	fjsA['mattr']['vorigin'] = 'somatic'
	fjsA['strandB'] = jsFile['strand']
	fjsA['strandA'] = jsFile['partner']['strand']
	if 'chimericreads' in jsFile:
		fjsA['clipreadB'] = jsFile['chimericreads'] if jsFile['chimericreads'] > 0 else 1
		if jsFile['ratio'] == 0:
 			 RatioB = 0.01
		else:
  			RatioB = jsFile['ratio']
		fjsA['totalreadB'] = int((fjsA['clipreadB']) // RatioB)
	if 'chimericreads' in jsFile['partner']:
		fjsA['clipreadA'] = jsFile['partner']['chimericreads'] if jsFile['partner']['chimericreads'] > 0 else 1
		if jsFile['partner']['ratio'] == 0:
			RatioA = 0.01
		else:
			RatioA = jsFile['partner']['ratio']
		fjsA['totalreadA'] = int((fjsA['clipreadA']) // RatioA)
	fjsB = fjsA.copy()
	posB = jsFile['position']
	posA = jsFile['partner']['position']
	fjsA['posA'] = posA
	fjsB['posB'] = posB
	fjsB['chrB'] = jsFile['chr']
	fjsA['chrA'] = jsFile['partner']['chr']
	if 'gene' in jsFile['partner']:
		fugeneA = jsFile['partner']['gene']
	else:
		fugeneA = fjsA['chrA']
	if 'gene' in jsFile:
		fugeneB = jsFile['gene']
	else:
		fugeneB = fjsB['chrB']
	fjsA = GETFUGENAB(fugeneA,fugeneB,fjsA)
	fjsB = GETFUGENAB(fugeneA,fugeneB,fjsB)
	#fjsA['fusiongene'] = fugeneA + ' > ' + fugeneB
	#fjsB['fusiongene'] = fugeneA + ' > ' + fugeneB
	return '\t'.join([fjsB['chrB'],str(posB),str(posB),json.dumps(fjsA,sort_keys=True)]) + '\n' + '\t'.join([fjsA['chrA'],str(posA),str(posA),json.dumps(fjsB,sort_keys=True)])+'\n'
def HANDLE_LIST(sam,pid,jsFile):
	FIN_O = ''
	for i in jsFile:
		fjsA = {}
		fjsA['dt'] = 2
		fjsA['sample'] = sam
		fjsA['mattr'] = {'project':'pcgp'}
		if pid:
			fjsA['mattr']['pmid'] = str(pid)
		fjsA['mattr']['vorigin'] = 'somatic'
		if 'type' in i:
			fjsA['svtype'] = i['type']
		fjsA['strandA'] = i['a']['strand']
		fjsA['strandB'] = i['b']['strand']
		if 'chimericreads' in i['a']:
			fjsA['clipreadA'] = i['a']['chimericreads'] if i['a']['chimericreads'] > 0 else 1
			if i['a']['ratio'] == 0:
				RatioA = 0.01
			else:
				RatioA = i['a']['ratio']
			fjsA['totalreadA'] = int((fjsA['clipreadA']) // RatioA)
		if 'chimericreads' in i['b']:
			fjsA['clipreadB'] = i['b']['chimericreads'] if i['b']['chimericreads'] > 0 else 1
			if i['b']['ratio'] == 0:
				RatioB = 0.01
			else:
				RatioB = i['b']['ratio']
			fjsA['totalreadB'] = int((fjsA['clipreadB']) // RatioB)
		fjsB = fjsA.copy()
		posA = i['a']['position']
		posB = i['b']['position']
		fjsA['posA'] = posA
		fjsB['posB'] = posB
		fjsA['chrA'] = i['a']['chr']
		fjsB['chrB'] = i['b']['chr']
		if 'name' in i['a']:
			fugeneA = i['a']['name']
		else:
			fugeneA = fjsA['chrA']
		if 'name' in i['b']:
			fugeneB = i['b']['name']
		else:
			fugeneB = fjsB['chrB']
		fjsA = GETFUGENAB(fugeneA,fugeneB,fjsA)
		fjsB = GETFUGENAB(fugeneA,fugeneB,fjsB)
		#fjsA['fusiongene'] = fugeneA + ' > ' + fugeneB
		#fjsB['fusiongene'] = fugeneA + ' > ' + fugeneB
		FIN_O += '\t'.join([fjsB['chrB'],str(posB),str(posB),json.dumps(fjsA,sort_keys=True)]) + '\n' + '\t'.join([fjsA['chrA'],str(posA),str(posA),json.dumps(fjsB,sort_keys=True)])+'\n'
	return FIN_O

for line in file:
	lineL = line.split('\t')
	js = json.loads(lineL[3])
	if 'codon' in line:
		out.write(HANDLE_DUX4IGH(lineL[0],lineL[4],js))
	elif 'position' not in line or lineL[6] == 'TARGET':
		continue
	elif 'typecode' in js:
		if js['typecode'] == 6:
			continue
		elif js['typecode'] == 9 and 'partner' in line:
			out.write(HANDLE_TY9(lineL[0],lineL[4],js))
		elif js['typecode'] == 8 and 'partner' in line:
			out.write(HANDLE_TY8(lineL[0],lineL[4],js))
	elif isinstance(js,list):
		out.write(HANDLE_LIST(lineL[0],lineL[4],js))
	else:
		print(line)
file.close()
out.close()
