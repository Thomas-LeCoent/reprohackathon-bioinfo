SRAID=Channel.from("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589", "SRR636531", "SRR636532", "SRR636533")

process SRA{
    publishDir "fastq/"
    input:
        val SRAID from SRAID
    output:
    file "${SRAID}_*.fastq.gz" into fastq
    script:
    """
    apt-get -y update
    apt-get -y install wget 
    wget -O ${SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${SRAID}/${SRAID}.1
    fastq-dump --gzip --split-files ${SRAID}.sra
    rm *.sra
    """
}
