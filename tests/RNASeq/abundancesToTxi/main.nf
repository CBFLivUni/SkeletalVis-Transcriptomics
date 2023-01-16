#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.species = "Human"
species="human"


include { abundancesToTxi } from '../../../modules/RNASeq/abundancesToTxi'

workflow test_abundancesToTxi {

  def input = []
    input = [file("tests/testData/test1.abundances"),file("tests/testData/test2.abundances")] 

    abundancesToTxi(input)
}
