/*
 * Include modules
 */

// SPAdes sreads
include { spades_sreads_assembly } from '../modules/ShortReads/spades_sreads.nf' params(outdir: params.outdir,
  threads: params.threads, spades_additional_parameters: params.spades_additional_parameters,
  shortreads_single: params.shortreads_single, shortreads_paired: params.shortreads_paired)

// Unicycler sreads
include { unicycler_sreads_assembly } from '../modules/ShortReads/unicycler_sreads.nf' params(outdir: params.outdir,
  threads: params.threads, unicycler_additional_parameters: params.unicycler_additional_parameters,
  shortreads_single: params.shortreads_single, shortreads_paired: params.shortreads_paired)

// Shovill sreads
include { shovill_sreads_assembly } from '../modules/ShortReads/shovill_sreads.nf' params(outdir: params.outdir,
  threads: params.threads, shovill_additional_parameters: params.shovill_additional_parameters,
  shortreads_single: params.shortreads_single, shortreads_paired: params.shortreads_paired)

/*
 * Module for assessing assembly qualities
 */
include { quast } from '../modules/QualityAssessment/quast.nf' params(threads: params.threads, outdir: params.outdir,
  longreads: params.longreads, shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)
include { multiqc } from '../modules/QualityAssessment/multiqc.nf' params(outdir: params.outdir)

workflow sreads_only_nf {
  take:
      preads
      sreads
  main:

  // Channels for quast
  spades_ch    = Channel.empty()
  unicycler_ch = Channel.empty()
  shovill_ch   = Channel.empty()

  // SPAdes
  if (!params.skip_spades) {
    spades_sreads_assembly(preads, sreads)
    spades_ch = spades_sreads_assembly.out[1]
  }
  // Unicycler
  if (!params.skip_unicycler) {
    unicycler_sreads_assembly(preads, sreads)
    unicycler_ch = unicycler_sreads_assembly.out[1]
  }
  // Shovill
  if (!params.skip_shovill && !params.shortreads_single && params.shortreads_paired) {
    shovill_sreads_assembly(preads)
    shovill_ch = shovill_sreads_assembly.out[1]
  }

  // Run quast
  quast(spades_ch.mix(unicycler_ch, shovill_ch), preads.concat(sreads).collect())

  // Run multiqc
  multiqc(quast.out[1].collect(), quast.out[2].distinct(), Channel.value('shortreads-only'))

}
