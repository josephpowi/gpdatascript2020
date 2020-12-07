if(process.argv.length!=4) {
	console.log('<cnv file> <sv file> output the same amount of cnv data to stdout')
	process.exit()
}

const cnvfile = process.argv[2]
const svfile = process.argv[3]

const fs=require('fs')


const sjid2name={}
{
	for(const line of fs.readFileSync('../../../files/hg19/nbl-hic/sjid2name',{encoding:'utf8'}).trim().split('\n')) {
		const l=line.split(' ')
		sjid2name[l[0]] = l[0]+'-'+l[1]
	}
}




const sample2sv = {}
/*
k: sample
   k: chr1
      v: [ { chr1/pos1/chr2/pos2/id } ]
*/

{
	const lines = fs.readFileSync(svfile,{encoding:'utf8'}).trim().split('\n')
/*
1	LiuYuTag	high
2	Sample(s)	SJNBL046422_C1
3	ChrA	1
4	PosA	2151235
5	OrientA	+
6	NumReadsA	6
7	ChrB	8
8	PosB	145148853
9	OrientB	-
10	NumReadsB	3
11	Type	CTX
33	CoverA	30
34	CoverB	42

*/
	for(let i=1; i<lines.length; i++) {
		const l = lines[i].split('\t')
		const sample =  sjid2name[l[1]]
		if(!sample2sv[sample]) sample2sv[sample]={}

		const chrA = 'chr'+l[3-1]
		const chrB = 'chr'+l[7-1]
		let v = Number.parseInt(l[4-1])
		if(Number.isNaN(v)) continue
		const posA=v-1
		v = Number.parseInt(l[8-1])
		if(Number.isNaN(v)) continue
		const posB=v-1
		const strandA=l[5-1]
		const strandB=l[9-1]
		const type=l[11-1]

		const cA=Number.parseInt(l[6-1])
		const cB=Number.parseInt(l[10-1])
		const tA=Number.parseInt(l[33-1])
		const tB=Number.parseInt(l[34-1])

		const j1={
			dt:5,
			sample:sample,
			strandA:strandA,
			chrB:chrB,
			posB:posB,
			strandB:strandB,
			clipreadA:cA,
			clipreadB:cB,
			totalreadA:cA,
			totalreadB:cB,
			svtype:type,
			mattr:{
				dna_assay:'wgs',
				vorigin:'somatic',
				project:'pedccl',
				pmid:'29284669'
			}
			}
		console.log(chrA+'\t'+posA+'\t'+posA+'\t'+JSON.stringify(j1))

		const j2={
			dt:5,
			sample:sample,
			strandB:strandB,
			chrA:chrA,
			posA:posA,
			strandA:strandA,
			clipreadA:cA,
			clipreadB:cB,
			totalreadA:cA,
			totalreadB:cB,
			svtype:type,
			mattr:{
				dna_assay:'wgs',
				vorigin:'somatic',
				project:'pedccl',
				pmid:'29284669'
			}
			}
		console.log(chrB+'\t'+posB+'\t'+posB+'\t'+JSON.stringify(j2))
	}
}




const lines=fs.readFileSync(cnvfile,{encoding:'utf8'}).trim().split('\n')

/*
1	sample	SJNBL046414_C1
2	chrom	chr1
3	loc.start	10101
4	loc.end	1860700
5	num.mark	9112
6	length.ratio	0.492
7	seg.mean	0.057

*/


for(let i=1; i<lines.length; i++) {
	const l = lines[i].split('\t')
	const v = Number.parseFloat(l[7-1])

	if(Number.isNaN(v)) {
		console.log('invalid logratio: '+lines[i])
		break
	}
	if(v==0) continue

	const sample = sjid2name[ l[0] ]

	const j={
		dt:4,
		sample:sample,
		value: v,
		mattr:{
			dna_assay:'wgs',
			vorigin:'somatic',
			project:'pedccl',
			pmid:'29284669'
		}
		}
	const chr = l[2-1]
	const start = Number.parseInt(l[3-1])-1
	const stop  = Number.parseInt(l[4-1])-1
	if(Number.isNaN(start) || Number.isNaN(stop)) continue
	console.log(chr+'\t'+start+'\t'+stop+'\t'+JSON.stringify(j))
}
