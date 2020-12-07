const fs=require('fs')


const files=[
	{sample:'SJNBL046414_C1',file:'/rgs01/resgen/prod/tartan/index/data/JinghuiZhang/CellBank/SJNBL046414_C1/WHOLE_GENOME/cnv-unpaired/SJNBL046414_C1_Single_CONSERTING_NoCREST_Mapability_100.txt'},
	{sample:'SJNBL046418_C1',file:'/rgs01/resgen/prod/tartan/index/data/JinghuiZhang/CellBank/SJNBL046418_C1/WHOLE_GENOME/cnv-unpaired/SJNBL046418_C1_Single_CONSERTING_NoCREST_Mapability_100.txt'},
	{sample:'SJNBL046420_C1',file:'/rgs01/resgen/prod/tartan/index/data/JinghuiZhang/CellBank/SJNBL046420_C1/WHOLE_GENOME/cnv-unpaired/SJNBL046420_C1_Single_CONSERTING_NoCREST_Mapability_100.txt'},
	{sample:'SJNBL046422_C1',file:'/rgs01/resgen/prod/tartan/index/data/JinghuiZhang/CellBank/SJNBL046422_C1/WHOLE_GENOME/cnv-unpaired/SJNBL046422_C1_Single_CONSERTING_NoCREST_Mapability_100.txt'},
	{sample:'SJNBL046424_C1',file:'/rgs01/resgen/prod/tartan/index/data/JinghuiZhang/CellBank/SJNBL046424_C1/WHOLE_GENOME/cnv-unpaired/SJNBL046424_C1_Single_CONSERTING_NoCREST_Mapability_100.txt'},
	{sample:'SJNBL046426_C1',file:'/rgs01/resgen/prod/tartan/index/data/JinghuiZhang/CellBank/SJNBL046426_C1/WHOLE_GENOME/cnv-unpaired/SJNBL046426_C1_Single_CONSERTING_NoCREST_Mapability_100.txt'}
	]


console.log('sample\tchrom\tloc.start\tloc.end\tnum.mark\tlength.ratio\tseg.mean')


for(const f of files) {
	const lines = fs.readFileSync(f.file,{encoding:'utf8'}).trim().split('\n')
	lines.shift()
	for(const s of lines) {
		console.log(f.sample+'\tchr'+s)
	}
}
