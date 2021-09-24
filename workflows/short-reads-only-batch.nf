/*
 * Include modules
 */

// SPAdes sreads
include { spades_batch } from '../modules/ShortReads/spades_sreads.nf'

// Unicycler sreads
// include { unicycler_batch } from '../modules/ShortReads/unicycler_sreads.nf'

// Shovill sreads
include { shovill_batch } from '../modules/ShortReads/shovill_sreads.nf'

/*
 * Module for assessing assembly qualities
 */
include { quast_batch } from '../modules/QualityAssessment/quast.nf'
include { multiqc_batch } from '../modules/QualityAssessment/multiqc.nf'

workflow sreads_only_batch_nf {

  take:
      input_tuple
  
  main:

  // Channels for quast
  spades_ch    = Channel.empty()
  unicycler_ch = Channel.empty()
  shovill_ch   = Channel.empty()

  // SPAdes
  if (!params.skip_spades) {
    spades_sreads_assembly(input_tuple)
    spades_ch = spades_sreads_assembly.out[1]
  }
  // Unicycler
  if (!params.skip_unicycler) {
    unicycler_sreads_assembly(input_tuple)
    unicycler_ch = unicycler_sreads_assembly.out[1]
  }
  // Shovill
  if (!params.skip_shovill) {
    shovill_batch(input_tuple.combine(Channel.from('spades', 'skesa', 'megahit')))
    shovill_ch = shovill_batch  .out[1]
  }

  // Get assemblies
  assemblies_ch = spades_ch.mix(unicycler_ch, shovill_ch)

  // Run quast
  quast_batch(assemblies_ch)

  // Run multiqc
  multiqc_batch(quast_batch.out[0].collect(), Channel.value("$workflow.runName"))

}
