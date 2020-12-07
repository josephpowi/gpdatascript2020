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


for(const line of fs.readFileSync('../../Pediatric/sampletable/target.samples',{encoding:'utf8'}).trim().split('\n')) {
/*
 * 1	diagnosis_group_short	HM
 * 2	diagnosis_group_full	Hematopoietic Malignancies
 * 3	diagnosis_short	BALL
 * 4	diagnosis_full	B-cell Acute Lymphoblastic Leukemia
 * 5	donor_name	SJCOGALL010887-PANKAK
 * 6	sample_name	SJCOGALL010887_D2-PANKAK
 * 7	sample_type	DIAGNOSIS
 * 8	target_id	PANKAK
 * 9	actual_RNA_SJID	SJCOGALL010887_D2
*/
	const l = line.split('\t')
	target2sample[l[8-1]+'-'+l[7-1]] = l[6-1]
}



const lines = fs.readFileSync('raw/LOH_every_case.txt',{encoding:'utf8'}).trim().split('\n')

/*
1	Sample	TARGET-10-CAAABC-09A-01D
2	chrom	chr1
3	loc.start	49342
4	loc.end	96091298
5	num.mark	41557
6	seg.mean	0.334
*/


const wrongsample=new Set()

for(let i=1; i<lines.length; i++) {

	const l = lines[i].split('\t')

	if(!targetbarcodeinuse.has(l[0])) {
		// excluded sample e.g. relapse
		wrongsample.add(l[0])
		continue
	}

	const targetid=l[0].split('_')[0] // TARGET-10-CAAABC-09A-01D


	const key = targetid.split('-')[2]+'-DIAGNOSIS'

	if(!target2sample[key]) {
		console.error('target LOH no sample for '+key)
		continue
	}

	const j={
		dt:10,
		sample:target2sample[key],
		segmean:Number.parseFloat(l[6-1]),
		mattr:{
			dna_assay:'cgi',
			project:'pantarget',
			vorigin:'somatic',
			pmid:'29489755',
		}
	}
	console.log(l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+JSON.stringify(j))
}


if(wrongsample.size) {
	console.error('target LOH wrong sample: '+[...wrongsample].join(' '))
}
