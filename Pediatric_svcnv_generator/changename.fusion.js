// some SV are missed by CGI (or no CGI available for that sample), add their fusions as SV





const fs=require('fs')
const exec=require('child_process').execSync










const target2sample = {}
const rnasjidmap = {}

for(const line of fs.readFileSync('../../Pediatric/sampletable/target.samples',{encoding:'utf8'}).trim().split('\n')) {
/*
1	diagnosis_group_short	ST
2	diagnosis_group_full	Solid Tumor
3	diagnosis_short	NBL
4	diagnosis_full	Neuroblastoma
5	donor_name	SJNBL017424-PASNPG
6	sample_name	SJNBL017424_D1-PASNPG
7	sample_type	DIAGNOSIS
8	target_id	PASNPG
9	actual_RNA_SJID	SJNBL017424_D1

*/
	const l = line.split('\t')
	target2sample[l[8-1]+'-'+l[7-1]] = l[6-1]

	if(l[9-1]) {
		rnasjidmap[l[9-1]]=l[6-1]
	}
}



const sample2fusion = new Map()
// key: sample
// val: {}
//    k: chra.posa.chrb.posb



{

	// no relpase samples in this file

	const lines = fs.readFileSync('../../../PanTarget/CGI_WXS_D/Fusion.txt',{encoding:'utf8'}).trim().split('\n')

	/*
	1	disease	AML
	2	sample	AML_diagnosis_PASWPD
	3	gene_a	UBTF
	4	chr_a	chr17
	5	position_a	42288315
	6	refseq_a	NM_014233
	7	strand_a	-
	8	gene_b	UBTF
	9	chr_b	chr17
	10	position_b	42289241
	11	refseq_b	NM_014233
	12	strand_b	-
	13	PLP_review	P
	14	PLP_geneA	UBTF
	15	PLP_geneB	
	16	source	added_05052017

	there could be duplicated records
	*/

	const dedup = new Set()

	for(let i=1; i<lines.length; i++) {

		const [disease, tarsample, geneA0, chrA, pstrA, ra, strandA, geneB0, chrB, pstrB, rb, strandB] = lines[i].split('\t')

		/*
		 * for OS_diagnosis_PAPFLB, remove all TP53 translocations except chr17-chr10
		 */
		if(tarsample=='OS_diagnosis_PAPFLB') {
			if(geneA0=='TP53' || geneB0=='TP53') {
				// is about TP53
				// following are two correct lines to keep for TP53, with highest score Liu Yu has said
				if(chrA=='chr17' && chrB=='chr10') {
				} else if(geneA0=='PMP22') {
				} else {
					// this is crapy line, jinghui doesn't want them to show up
					continue
				}
			}
		}

		const geneA = geneA0=='NA' ? null : geneA0
		const geneB = geneB0=='NA' ? null : geneB0


		const key = tarsample+chrA+pstrA+strandA+chrB+pstrB+strandB

		if(dedup.has(key)) continue
		dedup.add(key)

		const t = tarsample.split('_')
		const sampletype = t[1].toUpperCase() // must be upper case

		const sample= target2sample[ t[2]+'-'+sampletype ]

		if(!sample) {
			console.error('TARGET fusion sample not found: '+tarsample)
			continue
		}

		const posA=Number.parseInt(pstrA),
			posB=Number.parseInt(pstrB)

		if(Number.isNaN(posA) || Number.isNaN(posB)) {
			console.error('invalid coordinate: '+lines[i])
			continue
		}


		const fusiongene = (geneA || chrA) + ' > '+ (geneB || chrB)


		if(!sample2fusion.has(sample)) {
			sample2fusion.set(sample, new Set())
		}


		sample2fusion.get(sample).add( chrA+'.'+posA+'.'+strandA+'.'+chrB+'.'+posB+'.'+strandB )

		/*
		 * TALL RNA-seq are total, stranded
		 */
		const rnaseq_type = disease=='TALL' ? 'total' : 'polya'

		const j={
			dt:2,
			sample:sample,
			chrA:chrA,
			posA:posA,
			strandA:strandA,
			strandB:strandB,
			//fusiongene:fusiongene,
			mattr:{
				rna_assay: rnaseq_type,
				project:'pantarget',
				vorigin:'somatic',
				pmid:'29489755',
			}
		}
		if(geneA) j.geneA=geneA
		if(geneB) j.geneB=geneB

		console.log(chrB+'\t'+posB+'\t'+posB+'\t'+JSON.stringify(j))

		const k={
			dt:2,
			sample:sample,
			chrB:chrB,
			posB:posB,
			strandA:strandA,
			strandB:strandB,
			//fusiongene:fusiongene,
			mattr:{
				rna_assay: rnaseq_type,
				project:'pantarget',
				vorigin:'somatic',
				pmid:'29489755',
			}
		}
		if(geneA) k.geneA=geneA
		if(geneB) k.geneB=geneB
		console.log(chrA+'\t'+posA+'\t'+posA+'\t'+JSON.stringify(k))
	}
}




{
	const lines = fs.readFileSync('raw/NBL_final_fusions_manualReview_good.YLi',{encoding:'utf8'}).trim().split('\n')

	/*
	1	sample	SJNBL017066_D1
	2	YLi	bad
	3	Tag1	good
	4	Tag2	
	5	Tag3	show
	6	rating	good
	7	score	127.25
	8	geneA	ATP7A
	9	chrA	chrX
	10	posA	77227258
	11	ortA	+
	12	featureA	coding
	13	geneB	COX7B
	14	chrB	chrX
	15	posB	77154974
	16	ortB	+
	17	featureB	5utr
	*/

	let added=0

	for(let i=1; i<lines.length; i++) {

		const [ rnasjid, x1,x2,x3,x4,x5,x6, geneA0, chrA, pstrA, strandA, fa, geneB0, chrB, pstrB, strandB] = lines[i].split('\t')

		const geneA = geneA0=='NA' ? null : geneA0
		const geneB = geneB0=='NA' ? null : geneB0

		if(x1.startsWith('bad')) {
			// yongjin review
			continue
		}
		if(rnasjid.indexOf('_R')!=-1) {
			// no relapse from pan-target
			continue
		}

		const sample= rnasjidmap[rnasjid] || rnasjid

		const posA=Number.parseInt(pstrA),
			posB=Number.parseInt(pstrB)
		if(Number.isNaN(posA) || Number.isNaN(posB)) {
			console.error('NBL invalid coordinate: '+lines[i])
			continue
		}

		if( !sample2fusion.has( sample ) || !sample2fusion.get(sample).has( chrA+'.'+posA+'.'+strandA+'.'+chrB+'.'+posB+'.'+strandB)) {

			added++

			const fusiongene = (geneA || chrA) + ' > '+ (geneB || chrB)

			const j={
				dt:2,
				sample:sample,
				chrA:chrA,
				posA:posA,
				strandA:strandA,
				strandB:strandB,
				//fusiongene: fusiongene,
				mattr:{
					rna_assay:'polya',
					project:'pantarget',
					vorigin:'somatic',
					pmid:'29489755',
				}
			}
			if(geneA) j.geneA=geneA
			if(geneB) j.geneB=geneB
			console.log(chrB+'\t'+posB+'\t'+posB+'\t'+JSON.stringify(j))

			const k={
				dt:2,
				sample:sample,
				chrB:chrB,
				posB:posB,
				strandA:strandA,
				strandB:strandB,
				//fusiongene: fusiongene,
				mattr:{
					rna_assay:'polya',
					project:'pantarget',
					vorigin:'somatic',
					pmid:'29489755',
				}
			}
			if(geneA) k.geneA=geneA
			if(geneB) k.geneB=geneB
			console.log(chrA+'\t'+posA+'\t'+posA+'\t'+JSON.stringify(k))
		}
	}
	console.error('Added '+added+' target NBL fusion')
}





{
	const lines = fs.readFileSync('raw/WT_final_fusions_manualReview_goodones.YLi',{encoding:'utf8'}).trim().split('\n')

	/*
	1	sample	SJWLM019887_D1
	2	YLi	bad
	3	Tag1	good
	4	Tag2	
	5	Tag3	show
	6	rating	good
	7	score	148.26
	8	geneA	MAPKAPK5
	9	chrA	chr12
	10	posA	112308991
	11	ortA	+
	12	featureA	coding
	13	geneB	ACAD10
	14	chrB	chr12
	15	posB	112182445
	16	ortB	+
	17	featureB	coding
	18	sv_ort	>
	19	readsA	10
	20	readsB	15
	21	repeatA	0.069
	22	repeatB	0
	23	coverageA	157
	24	coverageB	490
	25	avg_coverage	237.81
	26	ratioA	0.044
	27	ratioB	0.1
	28	contig	TTTGCCAAGATTGACCAAGGTGACTTGATGACACCCCAGTTCACCCCTTATTATGTAGCACCCCAGGGCAAGCAAGCTCCACATATGCGGAACAAACTGGAAAGCTGACCGAATTTGTGTCTAACCTGGCGTGGGATTTCGCAGTCAAAGAAGGGTTCCGGGTTTTCAAAGAGATGCCCTTCACAAATCCGTTAACAAGGTCCTACCACACG
	29	type	INS
	30	frame	0,0,0,0,3,3
	31	sv_refseqA	NM_003668,NM_139078,NM_003668,NM_139078,,
	32	sv_refseqA_codon	,,,,,
	33	sv_refseqB	NM_001136538,NM_001136538,NM_025247,NM_025247,NM_001136538,NM_025247
	34	sv_refseqB_codon	603,603,572,572,603,000
	*/

	let added=0

	for(let i=1; i<lines.length; i++) {

		const [ rnasjid, x1,x2,x3,x4,x5,x6, geneA0, chrA, pstrA, strandA, fa, geneB0, chrB, pstrB, strandB ] = lines[i].split('\t')

		const geneA = geneA0=='NA' ? null : geneA0
		const geneB = geneB0=='NA' ? null : geneB0

		if(x1.startsWith('bad')) {
			// yongjin review
			continue
		}
		if(rnasjid.indexOf('_R')!=-1) {
			// no relapse from pan-target
			continue
		}

		const sample= rnasjidmap[rnasjid] || rnasjid

		const posA=Number.parseInt(pstrA)-1,
			posB=Number.parseInt(pstrB)-1
		if(Number.isNaN(posA) || Number.isNaN(posB)) {
			console.error('WT invalid coordinate: '+lines[i])
			continue
		}

		if( !sample2fusion.has( sample ) || !sample2fusion.get(sample).has( chrA+'.'+posA+'.'+strandA+'.'+chrB+'.'+posB+'.'+strandB)) {

			added++

			const fusiongene = (geneA || chrA) + ' > '+ (geneB || chrB)

			const j={
				dt:2,
				sample:sample,
				chrA:chrA,
				posA:posA,
				strandA:strandA,
				strandB:strandB,
				//fusiongene: fusiongene,
				mattr:{
					rna_assay:'polya',
					project:'pantarget',
					vorigin:'somatic',
					pmid:'29489755',
				}
			}
			if(geneA) j.geneA=geneA
			if(geneB) j.geneB=geneB
			console.log(chrB+'\t'+posB+'\t'+posB+'\t'+JSON.stringify(j))

			const k={
				dt:2,
				sample:sample,
				chrB:chrB,
				posB:posB,
				strandA:strandA,
				strandB:strandB,
				//fusiongene: fusiongene,
				mattr:{
					rna_assay:'polya',
					project:'pantarget',
					vorigin:'somatic',
					pmid:'29489755',
				}
			}
			if(geneA) k.geneA=geneA
			if(geneB) k.geneB=geneB
			console.log(chrA+'\t'+posA+'\t'+posA+'\t'+JSON.stringify(k))
		}
	}
	console.error('Added '+added+' target WT fusion')
}
