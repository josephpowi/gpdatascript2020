const fs=require('fs')

const lines=fs.readFileSync('raw/fusion_data_pre_2017-09-25_11_27_37.txt',{encoding:'utf8'}).trim().split('\n')
/*
11	Chr A	15
12	Pos A	81768973
13	Ori A	+
14	Locus A Soft Clips	7
15	Chr B	15
16	Pos B	81772165
17	Ori B	+
*/


console.log(lines[0])

const sample2data={}

for(let i=1; i<lines.length; i++) {
	const l=lines[i].split('\t')

	const sample=l[0]
	if(!sample2data[sample]) sample2data[sample]={}

	const k=l[11-1]+'.'+l[12-1]+'.'+l[13-1]+'.'+l[15-1]+'.'+l[16-1]+'.'+l[17-1]

	if(sample2data[sample][k]) continue

	sample2data[sample][k]=1

	console.log(lines[i])
}
