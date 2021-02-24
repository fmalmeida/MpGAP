/*
 * Include modules
 */

/*
 * Modules for assembling long reads
 */

// Canu assembler
include { canu_assembly } from '../modules/LongReads/canu.nf' params(outdir: params.outdir, lr_type: params.lr_type,
  canu_additional_parameters: params.canu_additional_parameters, threads: params.threads,
  genomeSize: params.genomeSize, longreads: params.longreads, shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single)

// Unicycler assembler
include { unicycler_lreads_assembly } from '../modules/LongReads/unicycler_lreads.nf' params(outdir: params.outdir,
  unicycler_additional_parameters: params.unicycler_additional_parameters, threads: params.threads, longreads: params.longreads,
  shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)

// Flye assembler
include { flye_assembly } from '../modules/LongReads/flye.nf' params(outdir: params.outdir, lr_type: params.lr_type,
  flye_additional_parameters: params.flye_additional_parameters, threads: params.threads, genomeSize: params.genomeSize,
  longreads: params.longreads, shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single)

// Raven assembler
include { raven_assembly } from '../modules/LongReads/raven.nf' params(outdir: params.outdir, threads: params.threads,
  raven_additional_parameters: params.raven_additional_parameters, longreads: params.longreads,
  shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)


/*
 * Modules for long reads assemblies polishment
 */

// Nanopolish (for nanopore data)
include { nanopolish } from '../modules/LongReads/nanopolish.nf' params(outdir: params.outdir, cpus: params.cpus, threads: params.threads,
  nanopolish_max_haplotypes: params.nanopolish_max_haplotypes, longreads: params.longreads, shortreads_paired: params.shortreads_paired,
  shortreads_single: params.shortreads_single, lr_type: params.lr_type)

// Medaka (for nanopore data)
include { medaka } from '../modules/LongReads/medaka.nf' params(medaka_sequencing_model: params.medaka_sequencing_model, threads: params.threads,
  outdir: params.outdir, longreads: params.longreads,
  shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)

// VariantCaller Pacbio
include { variantCaller } from '../modules/LongReads/variantCaller.nf' params(threads: params.threads, outdir: params.outdir, longreads: params.longreads,
  shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)

/*
 * Modules for Hybrid assemblies
 */

// Unicycler hybrid
include { unicycler_hybrid } from '../modules/Hybrid/unicycler_hybrid.nf' params(outdir: params.outdir,
  threads: params.threads, unicycler_additional_parameters: params.unicycler_additional_parameters,
  shortreads_single: params.shortreads_single, shortreads_paired: params.shortreads_paired)

// Unicycler hybrid
include { haslr_hybrid } from '../modules/Hybrid/haslr_hybrid.nf' params(outdir: params.outdir,
  threads: params.threads, haslr_additional_parameters: params.haslr_additional_parameters,
  shortreads_single: params.shortreads_single, shortreads_paired: params.shortreads_paired,
  longreads: params.longreads, lr_type: params.lr_type, genomeSize: params.genomeSize)

// SPAdes hybrid
include { spades_hybrid } from '../modules/Hybrid/spades_hybrid.nf' params(outdir: params.outdir,
  threads: params.threads, spades_additional_parameters: params.spades_additional_parameters,
  shortreads_single: params.shortreads_single, shortreads_paired: params.shortreads_paired,
  lr_type: params.lr_type)

// Pilon polish paired
include { pilon_polish } from '../modules/Hybrid/unicycler_polish.nf' params(outdir: params.outdir, threads: params.threads,
  pilon_memory_limit: params.pilon_memory_limit, shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single)

/*
 * Module for assessing assembly qualities
 */
include { quast } from '../modules/QualityAssessment/quast.nf' params(threads: params.threads, outdir: params.outdir, longreads: params.longreads,
  shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type, strategy_2: params.strategy_2)
include { multiqc } from '../modules/QualityAssessment/multiqc.nf' params(outdir: params.outdir)

workflow hybrid_nf {
  take:
      preads
      sreads
      lreads
      fast5
      fast5_dir
      bamFile
      nBams
  main:

      /*
       * Channels for placing the assemblies
       */
       // lreads
      lreads_assemblies_ch = Channel.empty()
      canu_ch              = Channel.empty()
      unicycler_ch         = Channel.empty()
      flye_ch              = Channel.empty()
      raven_ch             = Channel.empty()
       // hybrid
      hybrid_assemblies_ch = Channel.empty()
      spades_ch            = Channel.empty()
      unicycler_h_ch       = Channel.empty()
      haslr_ch             = Channel.empty()

      /*
       * Channels for placing polished assemblies
       */
      medaka_ch     = Channel.empty()
      nanopolish_ch = Channel.empty()
      arrow_ch      = Channel.empty()
      pilon_ch      = Channel.empty()

      /*
       * Full (default) hybrid mode
       */

      if (!params.strategy_2) {
        // SPAdes
        if (!params.skip_spades) {
          spades_hybrid(lreads, preads, sreads)
          spades_ch = spades_hybrid.out[1]
        }
        // Unicycler
        if (!params.skip_unicycler) {
          unicycler_hybrid(lreads, preads, sreads)
          unicycler_h_ch = unicycler_hybrid.out[1]
        }
        // Haslr
        if (!params.skip_haslr) {
          haslr_hybrid(lreads, preads, sreads)
          haslr_ch = haslr_hybrid.out[1]
        }

        // Get hybrid assemblies
        hybrid_assemblies_ch = spades_ch.mix(unicycler_h_ch, haslr_ch)
      }

      /*
       * Polish a long reads assembly
       */

      if (params.strategy_2) {
        /*
         * Canu
         */
        if (!params.skip_canu) {
          canu_assembly(lreads)
          canu_ch = canu_assembly.out[1]
        }

        /*
         * Flye
         */
        if (!params.skip_flye) {
          flye_assembly(lreads)
          flye_ch = flye_assembly.out[1]
        }

        /*
         * Unicycler
         */
        if (!params.skip_unicycler) {
          unicycler_lreads_assembly(lreads)
          unicycler_ch = unicycler_lreads_assembly.out[1]
        }

        /*
         * Raven
         */
        if (!params.skip_raven) {
          raven_assembly(lreads)
          raven_ch = raven_assembly.out[1]
        }

        // Get long reads assemblies
        lreads_assemblies_ch = canu_ch.mix(flye_ch, unicycler_ch, raven_ch)

        /*
         * Run medaka?
         */
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka(lreads_assemblies_ch.combine(lreads))
          medaka_ch = medaka.out[1]
        }

        /*
         * Run nanopolish?
         */
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish(lreads_assemblies_ch.combine(reads).combine(fast5).combine(fast5_dir))
          nanopolish_ch = nanopolish.out[0]
        }

        /*
         * VariantCaller?
         */
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller(lreads_assemblies_ch.combine(bamFile).combine(nBams))
          arrow_ch = variantCaller.out[1]
        }

        /*
         * Finally, run pilon for all
         */
        pilon_polish(
          lreads_assemblies_ch.mix(medaka_ch, nanopolish_ch, arrow_ch).combine(
            preads.combine(sreads).collect().toList()
          )
        )
        pilon_ch = pilon_polish.out[1]
      }

      // Run quast
      quast(
        hybrid_assemblies_ch.mix(lreads_assemblies_ch, pilon_ch).combine(
          preads.combine(sreads).collect().toList()
        )
      )

      // Run multiqc
      multiqc(quast.out[0].collect(), quast.out[1].distinct(), quast.out[2], Channel.value("$workflow.runName"))

}
