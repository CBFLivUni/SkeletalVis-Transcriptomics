process abundancesToTxi {

    tag "txi"
    publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

    label 'big_mem'
    input:
    path kallisto_abundance
    
    output:
    path "txi.RDS", emit: txiData
    path "txiScaledTPM.RDS"

    script:
    """
    Rscript ${params.scriptDir}/RNASeq/KallistoAbundancestoGeneCountMatrix.R --filepaths "$kallisto_abundance" --filenames "$kallisto_abundance" --species "${params.species}" --output "txi.RDS" --outputTPM "txiScaledTPM.RDS"

    """

}

