/* Test */

ALL_SRAID = "https://booster.pasteur.fr/static/files/primates/ref.nw.gz"

process downloadFastq{
    input: 
    val SRAID from ALL_SRAID

    // output: 
    // file "${SRAID}.sra" into fastq

    script:
    """
    echo ${SRAID}
    """
}   