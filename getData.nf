/* Test */

ID_LIST = Channel.from("SRR636567", "SRR636566")

process downloadFastq{
    input: 
    val SRAID from ID_LIST

    output: 
    file "${SRAID}.sra" into SRA_files

    script:
    """
    echo 'wget -O ${SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-7/${SRAID}/${SRAID}.1'
    """
}   
