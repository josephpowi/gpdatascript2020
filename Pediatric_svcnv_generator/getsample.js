if(process.argv.length!=3) {
	console.log('<in file>')
	process.exit()
}


const file=process.argv[2]

const fs=require('fs')
const readline=require('readline')

const r=readline.createInterface({input:fs.createReadStream(file,{encoding:'utf8'})})
const s=new Set()

r.on('line',line=>{
	s.add(JSON.parse(line.split('\t')[3]).sample)
})

r.on('close',()=>{
	console.log('#sample '+[...s].join(' '))
	console.error('GOT '+s.size+' samples')
})
