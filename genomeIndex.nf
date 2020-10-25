//Download the references .fasta for each human chromosome

ID=Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"MT")

process downloadHumanGenome{
    publishDir "fastagz/"
    input:
        val id from ID
    output:
        file "ref.fa" into fasta
    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${id}.fa.gz
    gunzip -c *.fa.gz >> ref.fa
    """
}

fasta.collectFile()

process createGenomeIndex{
    input:
        file ref.fa from fasta
        
    output:
        file "ref" into ref
        
    script:
    """
    mkdir ref
    STAR --runThreadN ${cpus} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ${ref.fa}
    """
}