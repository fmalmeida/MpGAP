/*
 * Modules for assembling long reads
 */

// Canu assembler
include { canu_assembly } from '../modules/LongReads/canu.nf' params(outdir: params.outdir, lr_type: params.lr_type,
  canu_additional_parameters: params.canu_additional_parameters, threads: params.threads,
  genomeSize: params.genomeSize)

// Unicycler assembler
include { unicycler_lreads } from '../modules/LongReads/unicycler_lreads.nf' params(outdir: params.outdir,
  unicycler_additional_parameters: params.unicycler_additional_parameters, threads: params.threads)

// Flye assembler
include { flye_assembly } from '../modules/LongReads/flye.nf' params(outdir: params.outdir, lr_type: params.lr_type,
  flye_additional_parameters: params.flye_additional_parameters, threads: params.threads,
  genomeSize: params.genomeSize)

// Raven assembler
include { raven_assembly } from '../modules/LongReads/raven.nf' params(outdir: params.outdir, threads: params.threads,
  raven_additional_parameters: params.raven_additional_parameters)

/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include { nanopolish } from '../modules/LongReads/nanopolish.nf' params(outdir: params.outdir, cpus: params.cpus, threads: params.threads,
                                                                                          nanopolish_max_haplotypes: params.nanopolish_max_haplotypes)

// Medaka (for nanopore data)
include { medaka } from '../modules/LongReads/medaka.nf' params(medaka_sequencing_model: params.medaka_sequencing_model, threads: params.threads, outdir: params.outdir)

// VariantCaller Pacbio
include { variantCaller } from '../modules/LongReads/variantCaller.nf' params(threads: params.threads, outdir: params.outdir)

/*
 * Module for assessing assembly qualities
 */
include { quast } from '../modules/QualityAssessment/quast.nf' params(threads: params.threads, outdir: params.outdir, longreads: params.longreads,
  shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)
include { multiqc } from '../modules/QualityAssessment/multiqc.nf' params(outdir: params.outdir)

workflow lreadsonly_nf {
  take:
      reads
      fast5
      fast5_dir
      bamFile
      nBams
  main:

      /*
       * Channels for placing the assemblies
       */
      canu_ch      = Channel.empty()
      unicycler_ch = Channel.empty()
      flye_ch      = Channel.empty()
      raven_ch     = Channel.empty()

      /*
       * Channels for placing polished assemblies
       */
      medaka_ch     = Channel.empty()
      nanopolish_ch = Channel.empty()
      arrow_ch      = Channel.empty()


      /*
       * Canu
       */
      if (!params.skip_canu) {
        canu_assembly(reads)
        canu_ch = canu_assembly.out[1]
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        flye_assembly(reads)
        flye_ch = flye_assembly.out[1]
      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        unicycler_lreads(reads)
        unicycler_ch = unicycler_lreads.out[1]
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        raven_assembly(reads)
        raven_ch = raven_assembly.out[1]
      }

      // gather assemblies
      assemblies_ch = canu_ch.mix(unicycler_ch, flye_ch, raven_ch)

      /*
       * Run medaka?
       */
      if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
        medaka(assemblies_ch.combine(reads))
        medaka_ch = medaka.out[1]
      }

      /*
       * Run nanopolish?
       */
      if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
        nanopolish(assemblies_ch.combine(reads, fast5, fast5_dir))
        nanopolish_ch = nanopolish.out[0]
      }

      /*
       * VariantCaller?
       */
      if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
        variantCaller(assemblies_ch.combine(bamFile, nBams))
        arrow_ch = variantCaller.out[1]
      }

      // Gather polishings
      polished_ch = medaka_ch.mix(nanopolish_ch, arrow_ch)

      /*
       * Run quast
       */
      quast(assemblies_ch.mix(polished_ch).combine(reads))

      /*
       * Run multiqc
       */
      multiqc(quast.out[1].collect(), quast.out[2].distinct(), Channel.value('longreads-only'), Channel.value("$workflow.runName"))
}
