/*
 * xiaotu's CGI CNV
 * TALL snp6 CNV
 *
 *
 * jian's WES CNV
 *
 *
 * record all sample names from CGI/TALL
 * do not show these samples from Jian's WES results
 *
 */



const fs=require('fs')


// record target barcodes in use
const targetbarcodeinuse = new Set()

/*
1	Cancer	BALL
2	sampletype	D
3	subtype	
4	Patient_ID	PAPDUV
5	CGI	1
6	WES	1
7	RNA	1
8	sjid	SJCOGALL010907_D4
9	wgs_d_tar	TARGET-10-PAPDUV-09A-01D
10	wxs_d_tar	TARGET-10-PAPDUV-09A-01D
11	wxs_d_sjid	SJCOGALL010907_D4
12	rna_d	TARGET-10-PAPDUV-09A-02R
13	rna_d_sjid	SJCOGALL010907_D3
*/
const lines1 = fs.readFileSync('../cohort/1699.DiagnosisOnly.pan_target_sample_info.update_Nov212016_OSbadAdded_update02082017_0A4I0S_back',{encoding:'utf8'}).trim().split('\n')
for(let i=1; i<lines1.length; i++) {
	const l = lines1[i].split('\t')
	const wgs_d_tar = l[9-1]
	if(wgs_d_tar) targetbarcodeinuse.add(wgs_d_tar)
	const wxs_d_tar = l[10-1]
	if(wxs_d_tar) targetbarcodeinuse.add(wxs_d_tar)
	const rna_d = l[12-1]
	if(rna_d) targetbarcodeinuse.add(rna_d)
}




const target2sample = {}
// k: PANKAK-diagnosis, v: SJCOGALL010887_D2-PANKAK


const validsamplenames = new Set()
// keep list of names like SJCOGALL010887_D2-PANKAK, so as to exclude relapse samples from WES CNV (which are not in the sample table)


for(const line of fs.readFileSync('../../Pediatric/sampletable/target.samples',{encoding:'utf8'}).trim().split('\n')) {
/*
1	diagnosis_group_short	HM
2	diagnosis_group_full	Hematopoietic Malignancies
3	diagnosis_short	BALL
4	diagnosis_full	B-cell Acute Lymphoblastic Leukemia
5	donor_name	SJCOGALL010887-PANKAK
6	sample_name	SJCOGALL010887_D2-PANKAK
7	sample_type	DIAGNOSIS
8	target_id	PANKAK
9	actual_RNA_SJID	SJCOGALL010887_D2
10	has_wgs	yes
11	has_wes	
12	has_rna	yes
13	target_DNA_barcode_diagnosis_wgs	TARGET-10-PANKAK-09A-01D
14	target_DNA_barcode_diagnosis_wes	
15	target_RNA_barcode_diagnosis	TARGET-10-PANKAK-09A-01R

*/
	const l = line.split('\t')
	const sample = l[6-1]

	target2sample[l[8-1]+'-'+l[7-1]] = sample

	validsamplenames.add(sample)
}





/*
 * record all sample names from CGI/TALL
 * do not show these samples from Jian's WES results
 */
const wesnotincludesamples = new Set()





const lines = fs.readFileSync('raw/CNV_every_case.txt',{encoding:'utf8'}).trim().split('\n')

/*
1	Sample	TARGET-10-CAAABC-09A-01D
2	chrchrom	chr1
3	loc.start	565401
4	loc.end	249209700
5	num.mark	216438
6	length.ratio	0.087
7	seg.mean	0.046
8	GMean	1.193
9	DMean	1.198
10	LogRatio	0.047
*/

const wrongsample=new Set()

for(let i=1; i<lines.length; i++) {

	const l = lines[i].split('\t')

	const targetid=l[0] // TARGET-10-CAAABC-09A-01D

	if(!targetbarcodeinuse.has(targetid)) {
		// should be an excluded sample, e.g. relapse
		wrongsample.add(targetid)
		continue
	}


	const key = targetid.split('-')[2]+'-DIAGNOSIS'

	const sample = target2sample[key]

	if(sample) {

		const j={
			dt:4,
			sample: sample,
			value:Number.parseFloat(l[10-1]),
			mattr:{
				dna_assay:'cgi',
				project:'pantarget',
				vorigin:'somatic',
				pmid:'29489755',
			}
		}
		console.log(l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+JSON.stringify(j))

		wesnotincludesamples.add( sample )

	} else {
		wrongsample.add( targetid )
		//console.error('target CNV no samplename for '+targetid)
	}
}

if(wrongsample.size) {
	console.error( 'target CNV: '+wrongsample.size+' samples wrong: '+[...wrongsample].join(' ') )
}





/* TALL snp6 CNV data
 * use USI to convert to sample name
 */

const usi2samplename = {}

/*
1	USI	PARASZ
16	Exome_id_D	SJALL015582_D1
29	sample_name	SJALL015582_D1-PARASZ

*/
for(const line of fs.readFileSync('../../Pediatric/sampletable/target.samples.tallsnp6array',{encoding:'utf8'}).trim().split('\n')) {
	const l = line.split('\t')
	usi2samplename[ l[0] ] = l[29-1]
}



const lines2 = fs.readFileSync('raw/TARGET_T_260pd_CBS_cnvFiltered_DGV_merged_hg19_IGV_w_genes.txt',{encoding:'utf8'}).trim().split('\n')

/*
1	TID	SJBALL001633_D1-TARGET-10-PASILW-09A-01D-snp6-vs-SJBALL001633_G1-TARGET-10-PASILW-10A-01D-snp6
2	SJID	SJBALL001633_D1_PASILW
3	UID	PASILW
4	chrom	4
5	cytoband	q25
6	loc.start	109041477
7	loc.end	109060364
8	num.mark	18
9	seg.mean	-0.7375
*/

for(let i=1; i<lines2.length; i++) {

	const l = lines2[i].split('\t')

	const chr=l[4-1]
	if(!chr) continue

	const t=l[1].split('_')

	const usi = l[3-1]

	const  samplename = usi2samplename[usi]
	if(!samplename) {
		console.error('unknown usi from TALL snp6 cnv: '+usi)
		continue
	}

	const j={
		dt:4,
		sample:samplename,
		value:Number.parseFloat(l[9-1]),
		mattr:{
			dna_assay:'snp6',
			project:'pantarget',
			vorigin:'somatic',
			pmid:'28671688',
		}
	}

	console.log('chr'+l[4-1]+'\t'+l[6-1]+'\t'+l[7-1]+'\t'+JSON.stringify(j))

	wesnotincludesamples.add( samplename )
}





// WES CNV from jian
// diagnosis tumor only, no relapse from this file
//
//

const skippedsamples = new Set()
const wesusesamples = new Set()

for(const line of fs.readFileSync('raw/wes.conserting.cnv', {encoding:'utf8'}).trim().split('\n')) {

	const l = line.split('\t')
	const j = JSON.parse( l[3] )

	if(wesnotincludesamples.has(j.sample)) {
		skippedsamples.add( j.sample )
		continue
	}

	if(!validsamplenames.has(j.sample)) {
		skippedsamples.add(j.sample)
		continue
	}

	wesusesamples.add(j.sample)
	j.mattr = {
		dna_assay:'wes',
		project:'pantarget',
		vorigin:'somatic',
		pmid:'29489755',
	}
	console.log( l[0]+'\t'+l[1]+'\t'+l[2]+'\t'+JSON.stringify( j ) )
}
console.error( 'WES CNV use '+wesusesamples.size+' samples, skipped '+skippedsamples.size )
