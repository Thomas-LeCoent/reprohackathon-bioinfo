//accesion id
//SRAID=Channel.from("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589", "SRR636531", "SRR636532", "SRR636533")
SRAID=Channel.from("SRR628585")

//Download SRA files
process downloadSRA{
	publishDir "sraFiles/"
	input:
	val id from SRAID
	output:
	file "${id}.sra" into sraF

	script:
	"""
	wget -O ${id}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${id}/${id}.1
	"""
}

//Get fastq files from the sra
process fastqDump{

    publishDir "fastq/"
    
    input:
	file sraFile from sraF
	
    output:
	file "*.fastq.gz" into fastq

    script:
    """
    fasterq-dump --split-files ${sraFile}
    gzip *.fastq
    """
}

//Download annotation in .gtf
process downloadGenomeAnnotation{

    publishDir "gtf/"
    
    output:
    file "annot.gtf" into gtf
    
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

    publishDir "fastagz/"
    
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
    
cpus = 16

process createGenomeIndex{

	publishDir "genomeIndex/"
	
	input:
	file (genome) from fasta.collectFile() //Put all the extracted files in a single one
	
	output:
	path "ref" into index_chan 

	script:
	"""
	STAR --runThreadN ${cpus} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ${genome}
	"""
}

fastq = Channel.fromFilePairs('fastq/*_{1,2}.fastq.gz', flat:true)
fastq_2=fast_q.combine(index_chan)
process mapping{
	publishDir "BAM/"
	
	input:
	tuple val(sampleID), file(r1), file(r2), path(r3) from fastq_2

	output:
	file "${sampleID}.bam" into bam_chan

	script:
	"""
	STAR --outSAMstrandField intronMotif \
	--outFilterMismatchNmax 4 \
	--outFilterMultimapNmax 10 \
	--genomeDir ${r3} \
	--readFilesIn <(gunzip -c ${r1}) <(gunzip -c ${r2}) \
	--runThreadN 16 \
	--outSAMunmapped None \
	--outSAMtype BAM SortedByCoordinate \
	--outStd BAM_SortedByCoordinate \
	--genomeLoad NoSharedMemory \
	--limitBAMsortRAM 100000000 \
	> ${sampleID}.bam
	"""
}
process samtools{
	publishDir "bam/"
	
	input:
	file(bam_to_index) from bam_chan

	output:
	file "${bam_to_index}.bai" into end
	file "$bam_to_index" into ready_to_count
	
	script:
	"""
	samtools index $bam_to_index > ${bam_to_index.baseName}.bai
	"""
}

process read_count{
	publishDir "read_count/"

	input:
	file(bam) from ready_to_count
	file(gtf) from gtf

	output:
	file "${bam.baseName}.counts" into read_count

	script:
	"""
	featureCounts $bam -T $cpus -t gene -g gene_id -s 0 -a $gtf -o ${bam.baseName}.counts
	"""
}
