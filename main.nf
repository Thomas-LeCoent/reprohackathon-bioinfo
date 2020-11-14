//accesion id
SRAID = Channel.from("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
// SRAID=Channel.from("SRR628585")

//Download SRA files
process downloadSRA{
	publishDir "files/downloadSRA/"

	input:
	val sraid from SRAID

	output:
		tuple val(sraid), file("${sraid}.sra") into sraFiles

	script:
	"""
	wget -O ${sraid}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${sraid}/${sraid}.1
	"""
}

// Get fastq files from the sra
process fastqDump{
    publishDir "files/fastqDump/"
    
    input:
		tuple val(sraid), file(sra) from sraFiles
	
    output:
		tuple val(sraid), file("*_1.fastq.gz"),file("*_2.fastq.gz") into fastq
	
    script:
    """
    fasterq-dump --split-files ${sra}
    gzip *.fastq
    """
}


//Download annotation in .gtf
process downloadGenomeAnnotation{

    publishDir "files/downloadGenomeAnnotation/"
    
    //input : val gtfURL
    
    output:
    file "annot.gtf" into gtf
    //
    
    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip -c *.gtf.gz > annot.gtf
    """
}

//human chromosomes
ID=Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"MT")

//Download the references .fasta
process downloadHumanGenome{

    publishDir "files/downloadHumanGenome/"
    
    input:
        val id from ID
	
    output:
        file "ref.fa" into fasta
	
    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${id}.fa.gz
    gunzip -c *fa.gz > ref.fa
    """
}
    

process createGenomeIndex{

	publishDir "files/createGenomeIndex/"
	
	input:
	file (genome) from fasta.collectFile()
	
	output:
	path "ref" into index_chan 

	script:
	"""
	STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ${genome}
	"""
}

//fastq = Channel.fromFilePairs('fastq/*_{1,2}.fastq.gz', flat:true)
//fastq=fastq.combine(index_chan)//combine index and fastq for the mapping
process mapping{
	publishDir "files/mapping/"
	
	input:
	tuple val(sraid), file(r1), file(r2) from fastq
	// path index from index_chan
	path index from index_chan.first()

	output:
	file "${sraid}.bam" into bam_chan, bam_chan2

	script:
	"""
	STAR --outSAMstrandField intronMotif \
	--outFilterMismatchNmax 4 \
	--outFilterMultimapNmax 10 \
	--genomeDir ${index} \
	--readFilesIn <(gunzip -c ${r1}) <(gunzip -c ${r2}) \
	--runThreadN ${task.cpus} \
	--outSAMunmapped None \
	--outSAMtype BAM SortedByCoordinate \
	--outStd BAM_SortedByCoordinate \
	--genomeLoad NoSharedMemory \
	--limitBAMsortRAM ${task.memory.toBytes()} \
	> ${sraid}.bam
	"""
}

process samtools{
	publishDir "files/samtools/"
	
	input:
	file(bam_to_index) from bam_chan

	output:
	file "${bam_to_index}.bai" into end
	
	
	script:
	"""
	samtools index $bam_to_index 
	"""
}

process featureCounts{
	publishDir "files/featureCounts/"

	input:
	file(bam) from bam_chan2.collect()
	file(gtf) from gtf

	output:
	file "counts.txt" into read_count

	script:
	"""
	featureCounts $bam -T ${task.cpus} -t gene -g gene_id -s 0 -a $gtf -o counts.txt
	"""
}


/* process DESeq{

	input:
	file count from read_count
	
	output:
	file "Result_DESeq.txt"
	
	script:
	template "DESeq2_count.R"

}
*/

workflow.onComplete = {
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}


workflow.onError = {
    println "Oops .. something when wrong"
}
