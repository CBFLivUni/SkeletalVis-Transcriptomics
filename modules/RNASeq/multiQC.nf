process multiQC {

    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    label 'standard'
    input:
    path ('fastqc/*')
    path ('kallisto/*')
    path ('fastqscreen/*')

    output:
    path "*multiqc_report.html"
    path "*_data"

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

    multiqc .
    """

}

