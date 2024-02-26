#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = "/home/diagnostics/CodecRun/codecPipeline/sample_sheet.csv"
process splitLetters {
	output: 
	stdout
	"""
	cat  ${params.input}
	"""
}


workflow {
	Channel
		
		.fromPath(params.input) 
		.splitCsv(header:false)
		.flatten() 
		.view()
		.map { it}
		.set {Sample}
}
