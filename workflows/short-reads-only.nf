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
include { quast as quast_sreads_spades; quast as quast_sreads_unicycler ; quast as quast_sreads_shovill } \
       from '../modules/QualityAssessment/quast.nf' params(threads: params.threads, outdir: params.outdir, longreads: params.longreads,
            shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)

workflow sreads_only_nf {
  take:
      preads
      sreads
  main:

  // SPAdes
  if (!params.skip_spades) {
    spades_sreads_assembly(preads, sreads)
    quast_sreads_spades(spades_sreads_assembly.out[1], preads.concat(sreads).collect())
  }
  // Unicycler
  if (!params.skip_unicycler) {
    unicycler_sreads_assembly(preads, sreads)
    quast_sreads_unicycler(unicycler_sreads_assembly.out[1], preads.concat(sreads).collect())
  }
  // Shovill
  if (!params.skip_shovill && !params.shortreads_single && params.shortreads_paired) {
    shovill_sreads_assembly(preads)
    quast_sreads_shovill(shovill_sreads_assembly.out[1], preads.concat(sreads).collect())
  }

}
