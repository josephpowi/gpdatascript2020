const fs=require('fs')


const target2sample = {}


for(const line of fs.readFileSync('../../Pediatric/sampletable/target.samples',{encoding:'utf8'}).trim().split('\n')) {
/*
1	diagnosis_group_short	HM
2	diagnosis_group_full	Hematopoietic Malignancies
3	diagnosis_short	TALL
4	diagnosis_full	T-cell Acute Lymphoblastic Leukemia
5	donor_name	CAAABC
6	sample_name	CAAABC-diagnosis
7	sample_type	DIAGNOSIS
8	target_id	CAAABC

*/
	const l = line.split('\t')
	target2sample[l[8-1]+'-'+l[7-1]] = l[6-1]
}




// this file has only diagnosis samples

const lines = fs.readFileSync('raw/SV_table.txt',{encoding:'utf8'}).trim().split('\n')

/*
 
### old file, uncleaned, with relapse samples
2	D_R	D
4	FusionBuilder_gene	CHST13,GRK1
29	>Id	TARGET-40-PARKAF-01A-01D
30	LeftChr	chr1
31	LeftPosition	118121700
32	LeftStrand	-
33	LeftLength	534
34	RightChr	chr4
35	RightPosition	164697661
36	RightStrand	-
55	Type	interchromosomal






1	cancer	OS
2	patient	0A4I9K
3	sampletype	diagnosis
4	gene_a	
5	refseq_a	
6	chr_a	1
7	position_a	118331947
8	strand_a	-
9	gene_b	SNX27
10	refseq_b	NM_030918
11	chr_b	1
12	position_b	151647456
13	strand_b	+
14	Discordant_readcount	9
15	PLP_review	B
16	PLP_geneA	
17	PLP_geneB	

*/


let grk1lines=0
let nbaslines_nbl=0
let nbaslines_wt=0


for(let i=1; i<lines.length; i++) {

	const [cancer,patient,sampletype,genea0,refseqa,chra0,posastr,stranda,geneb0,refseqb,chrb0,posbstr,strandb,discordantreadstr] = lines[i].split('\t')

	const chra = 'chr'+chra0
	const chrb = 'chr'+chrb0

	// genea and geneb names may be replaced
	let genea = genea0
	let geneb = geneb0

	/*
	drop GRK1 - CHST13
	*/
	if((genea=='GRK1' && geneb=='CHST13') || (genea=='CHST13' && geneb=='GRK1')) {
		grk1lines++
		continue
	}

	/*
	 * for nbl & wt, replace NBAS with MYCN
	if(cancer=='NBL') {

		let replaced=false
		if(genea=='NBAS') {
			genea='MYCN'
			replaced=true
		}
		if(geneb=='NBAS') {
			geneb='MYCN'
			replaced=true
		}
		if(replaced) nbaslines_nbl++

	} else if(cancer=='WT') {

		let replaced=false
		if(genea=='NBAS') {
			genea='MYCN'
			replaced=true
		}
		if(geneb=='NBAS') {
			geneb='MYCN'
			replaced=true
		}
		if(replaced) nbaslines_wt++

	}

	 */





	const target = patient+'-'+sampletype.toUpperCase()

	const samplename = target2sample[ target ]


	if(!samplename) {
		console.error('target SV wrong sample: '+target)
		continue
	}

	const posa=Number.parseInt(posastr)
	const posb=Number.parseInt(posbstr)

/*
	let type
	switch(l[55-1]) {
	case 'artifact':
		type='artifact'
		break
	case 'complex':
		type='complex'
		break
	case 'deletion':
		type='DEL'
		break
	case 'distal-duplication':
	case 'distal-duplication-by-mobile-element':
	case 'tandem-duplication':
		type='INS'
		break
	case 'interchromosomal':
		type='CTX'
		break
	case 'inversion':
	case 'probable-inversion':
		type='INV'
		break
	default:
		console.error('ERR: unknown target SV type: '+l[55-1])
	}
	*/




	const j1={
		dt:5,
		sample: samplename,
		strandA:stranda,
		chrB:chrb,
		posB:posb,
		strandB:strandb,
		//svtype:type,
		mattr: {
			dna_assay:'cgi',
			project:'pantarget',
			vorigin:'somatic',
			pmid:'29489755',
		}
	}
	if(genea!='') j1.geneA = genea
	if(geneb!='') j1.geneB = geneb
	console.log(chra+'\t'+posa+'\t'+posa+'\t'+JSON.stringify(j1))



	const j2={
		dt:5,
		sample: samplename,
		strandB:strandb,
		chrA:chra,
		posA:posa,
		strandA:stranda,
		//svtype:type,
		mattr:{
			dna_assay:'cgi',
			project:'pantarget',
			vorigin:'somatic',
			pmid:'29489755',
		}
	}
	if(genea!='') j2.geneA = genea
	if(geneb!='') j2.geneB = geneb
	console.log(chrb+'\t'+posb+'\t'+posb+'\t'+JSON.stringify(j2))
}

console.error('pan-TARGET SV: dropped '+grk1lines+' lines of GRK1-CHST13')
console.error('pan-TARGET SV: '+nbaslines_nbl+' lines of NBAS replaced in NBL')
console.error('pan-TARGET SV: '+nbaslines_wt+' lines of NBAS replaced in WT')
