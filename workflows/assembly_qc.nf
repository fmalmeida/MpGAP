/*
 * Include modules
 */

/*
 * Module for assessing assembly qualities
 */
include { quast   } from '../modules/local/QualityAssessment/quast.nf'
include { multiqc } from '../modules/local/QualityAssessment/multiqc.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'  

workflow ASSEMBLY_QC {
  take:
    input_tuple
    ch_versions
  
  main:

    // Run quast (with all)
    quast(input_tuple)
    ch_versions = ch_versions.mix(quast.out.versions.first())

    // collect software versions
    CUSTOM_DUMPSOFTWAREVERSIONS( ch_versions.unique().collectFile(name: 'collated_versions.yml') )

    // Run multiqc
    multiqc(
      quast.out.results.groupTuple(by: [0,1,2]),
      Channel.value("$workflow.runName"),
      CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect().first()
    )

}
