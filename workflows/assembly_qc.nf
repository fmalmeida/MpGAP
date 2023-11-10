/*
 * Include modules
 */

/*
 * Module for assessing assembly qualities
 */
include { quast   } from '../modules/QualityAssessment/quast.nf'
include { multiqc } from '../modules/QualityAssessment/multiqc.nf'

workflow ASSEMBLY_QC {
  take:
    input_tuple
  
  main:

    // Run quast (with all)
    quast(input_tuple)

    // Run multiqc
    multiqc(quast.out.results.groupTuple(by: [0,1,2]), Channel.value("$workflow.runName"))

}
