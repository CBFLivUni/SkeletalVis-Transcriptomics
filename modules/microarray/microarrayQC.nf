process microarrayQC {

    tag "qc"
    publishDir path: "${params.outdir}/qc", mode: 'copy'

    label 'short_job'

    input:
    val qcData

    output:
    path "*"

    when:
    params.platform != "2C-Agilent"

    script:
    
    """
    Rscript ${params.scriptDir}/microarray/microarrayQC.R --inputfile $qcData --comparisonsTable "${params.comparisonsTable}" --platform "${params.platform}" --accessionNumber "${params.accessionNumber}" --offset "${params.offset}" --outputFile "${params.accessionNumber}"_qc.html --outputDirectory ./
    """

}

