//Download the references .fasta for each human chromosome

ChrNames = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"MT")
cpus = 8

process downloadHumanChromosomes{
    input:
        val chr from ChrNames
    output:
        file "Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa" into ChrFiles
    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz
    gunzip *.fa.gz > ${chr}.fa
    """
}


ChrFiles.collectFile().println()

/*
process createGenomeIndex{
    container="evolbioinfo/star:v2.7.6a"
    
    input:
        file (genome) from fasta.collectFile()
        
    script:
    """
    mkdir ref
    STAR --runThreadN ${cpus} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ${genome}
    """
}*/
