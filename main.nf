//param
description=Channel.fromPath("./descriptionMutation.txt")

//accesion id
SRAID = Channel.from("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")




process Deseq2{
	publishDir "fileDESeq/Results/"
	input:
	file count from read_count
	file des from description
	
	output:
	file("results.txt") into fin
	
	script:
	template "DESeq2_count.R"
}


workflow.onComplete = {
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}


workflow.onError = {
    println "Oops .. something when wrong"
}
