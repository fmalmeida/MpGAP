/*
 * Include modules
 */

/*
 * Module for prefix evaluation
 */
include { define_prefix } from '../modules/misc/define_prefix.nf'

// SPAdes sreads
include { spades_sreads_assembly } from '../modules/ShortReads/spades_sreads.nf'

// Unicycler sreads
include { unicycler_sreads_assembly } from '../modules/ShortReads/unicycler_sreads.nf'

// Shovill sreads
include { shovill_sreads_assembly } from '../modules/ShortReads/shovill_sreads.nf'

/*
 * Module for assessing assembly qualities
 */
include { quast } from '../modules/QualityAssessment/quast.nf'
include { multiqc } from '../modules/QualityAssessment/multiqc.nf'

workflow sreads_only_nf {
  take:
      preads
      sreads
  main:

  /*
   * Check input reads in order to evaluate better prefix
   */
  define_prefix(Channel.value(''), preads, sreads)
  prefix_ch = define_prefix.out[0]

  // Channels for quast
  spades_ch    = Channel.empty()
  unicycler_ch = Channel.empty()
  shovill_ch   = Channel.empty()

  // SPAdes
  if (!params.skip_spades) {
    spades_sreads_assembly(preads.combine(sreads).combine(prefix_ch))
    spades_ch = spades_sreads_assembly.out[1]
  }
  // Unicycler
  if (!params.skip_unicycler) {
    unicycler_sreads_assembly(preads.combine(sreads).combine(prefix_ch))
    unicycler_ch = unicycler_sreads_assembly.out[1]
  }
  // Shovill
  if (!params.skip_shovill && !params.shortreads_single && params.shortreads_paired) {
    preads.combine(prefix_ch).combine(Channel.value(['spades', 'skesa', 'megahit'])).view(0)
    shovill_sreads_assembly(preads.combine(prefix_ch).combine(Channel.value(['spades', 'skesa', 'megahit'])))
    shovill_ch = shovill_sreads_assembly.out[1]
  }

  // Get assemblies
  assemblies_ch = spades_ch.mix(unicycler_ch, shovill_ch)

  // Run quast
  quast(
    assemblies_ch.combine(preads.combine(sreads).collect().toList()).combine(prefix_ch)
  )

  // Run multiqc
  multiqc(quast.out[0].collect(), prefix_ch, Channel.value("$workflow.runName"))

}
