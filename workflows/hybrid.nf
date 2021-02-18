/*
 * Include modules
 */

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
include { nanopolish as nanopolish_canu;
          nanopolish as nanopolish_unicycler;
          nanopolish as nanopolish_flye;
          nanopolish as nanopolish_raven } from '../modules/LongReads/nanopolish.nf' params(outdir: params.outdir, cpus: params.cpus, threads: params.threads,
                                                                                          nanopolish_max_haplotypes: params.nanopolish_max_haplotypes)

// Medaka (for nanopore data)
include { medaka as medaka_canu;
          medaka as medaka_unicycler;
          medaka as medaka_flye;
          medaka as medaka_raven } from '../modules/LongReads/medaka.nf' params(medaka_sequencing_model: params.medaka_sequencing_model, threads: params.threads, outdir: params.outdir)

// VariantCaller Pacbio
include { variantCaller as variantCaller_canu;
          variantCaller as variantCaller_unicycler;
          variantCaller as variantCaller_flye;
          variantCaller as variantCaller_raven } from '../modules/LongReads/variantCaller.nf' params(threads: params.threads, outdir: params.outdir)

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
include { pilon_polish as pilon_polish_flye; pilon_polish as pilon_polish_flye_nanopolish;
          pilon_polish as pilon_polish_flye_medaka; pilon_polish as pilon_polish_flye_variantCaller;
          pilon_polish as pilon_polish_raven; pilon_polish as pilon_polish_raven_nanopolish;
          pilon_polish as pilon_polish_raven_medaka; pilon_polish as pilon_polish_raven_variantCaller;
          pilon_polish as pilon_polish_canu; pilon_polish as pilon_polish_canu_nanopolish;
          pilon_polish as pilon_polish_canu_medaka; pilon_polish as pilon_polish_canu_variantCaller;
          pilon_polish as pilon_polish_unicycler; pilon_polish as pilon_polish_unicycler_nanopolish;
          pilon_polish as pilon_polish_unicycler_medaka; pilon_polish as pilon_polish_unicycler_variantCaller } \
\
          from '../modules/Hybrid/unicycler_polish.nf' params(outdir: params.outdir, threads: params.threads,
            pilon_memory_limit: params.pilon_memory_limit, shortreads_paired: params.shortreads_paired,
            shortreads_single: params.shortreads_single)

/*
 * Module for assessing assembly qualities
 */
include { quast } from '../modules/QualityAssessment/quast.nf' params(threads: params.threads, outdir: params.outdir, longreads: params.longreads,
  shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)
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
       * Quast channels
       */

      // Hybrid st 1
      spades_ch = Channel.empty()
      unicycler_hybrid_ch = Channel.empty()
      haslr_ch        = Channel.empty()

      // Hybrid st 2
      canu_ch         = Channel.empty()
      nanopol_canu_ch = Channel.empty()
      medaka_canu_ch  = Channel.empty()
      arrow_canu_ch   = Channel.empty()
      raven_ch         = Channel.empty()
      nanopol_raven_ch = Channel.empty()
      medaka_raven_ch  = Channel.empty()
      arrow_raven_ch   = Channel.empty()
      flye_ch         = Channel.empty()
      nanopol_flye_ch = Channel.empty()
      medaka_flye_ch  = Channel.empty()
      arrow_flye_ch   = Channel.empty()
      unicycler_ch         = Channel.empty()
      nanopol_unicycler_ch = Channel.empty()
      medaka_unicycler_ch  = Channel.empty()
      arrow_unicycler_ch   = Channel.empty()

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
          unicycler_hybrid_ch = unicycler_hybrid.out[1]
        }
        // Haslr
        if (!params.skip_haslr) {
          haslr_hybrid(lreads, preads, sreads)
          haslr_ch = haslr_hybrid.out[1]
        }
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
          pilon_polish_canu(canu_assembly.out[1], preads.concat(sreads).collect())
          canu_ch = pilon_polish_canu.out[1]

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_canu(canu_assembly.out[1], lreads, fast5, fast5_dir)
            pilon_polish_canu_nanopolish(nanopolish_canu.out[0], preads.concat(sreads).collect())
            nanopol_canu_ch = pilon_polish_canu_nanopolish.out[1]
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_canu(canu_assembly.out[1], lreads)
            pilon_polish_canu_medaka(medaka_canu.out[1], preads.concat(sreads).collect())
            medaka_canu_ch = pilon_polish_canu_medaka.out[1]
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_canu(canu_assembly.out[1], bamFile, nBams)
            pilon_polish_canu_variantCaller(variantCaller_canu.out[1], preads.concat(sreads).collect())
            arrow_canu_ch = pilon_polish_canu_variantCaller.out[1]
          }
        }

        /*
         * Flye
         */
        if (!params.skip_flye) {
          flye_assembly(lreads)
          pilon_polish_flye(flye_assembly.out[1], preads.concat(sreads).collect())
          flye_ch = pilon_polish_flye.out[1]

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_flye(flye_assembly.out[1], lreads, fast5, fast5_dir)
            pilon_polish_flye_nanopolish(nanopolish_flye.out[0], preads.concat(sreads).collect())
            nanopol_flye_ch = pilon_polish_flye_nanopolish.out[1]
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_flye(flye_assembly.out[1], lreads)
            pilon_polish_flye_medaka(medaka_flye.out[1], preads.concat(sreads).collect())
            medaka_flye_ch = pilon_polish_flye_medaka.out[1]
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_flye(flye_assembly.out[1], bamFile, nBams)
            pilon_polish_flye_variantCaller(variantCaller_flye.out[1], preads.concat(sreads).collect())
            arrow_flye_ch = pilon_polish_flye_variantCaller.out[1]
          }

        }

        /*
         * Unicycler
         */
        if (!params.skip_unicycler) {
          unicycler_lreads(lreads)
          pilon_polish_unicycler(unicycler_lreads.out[1], preads.concat(sreads).collect())
          unicycler_ch = pilon_polish_unicycler.out[1]

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_unicycler(unicycler_lreads.out[1], lreads, fast5, fast5_dir)
            pilon_polish_unicycler_nanopolish(nanopolish_unicycler.out[0], preads.concat(sreads).collect())
            nanopol_unicycler_ch = pilon_polish_unicycler_nanopolish.out[1]
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_unicycler(unicycler_lreads.out[1], lreads)
            pilon_polish_unicycler_medaka(medaka_unicycler.out[1], preads.concat(sreads).collect())
            medaka_unicycler_ch = pilon_polish_unicycler_medaka.out[1]
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_unicycler(unicycler_lreads.out[1], bamFile, nBams)
            pilon_polish_unicycler_variantCaller(variantCaller_unicycler.out[1], preads.concat(sreads).collect())
            arrow_unicycler_ch = pilon_polish_unicycler_variantCaller.out[1]
          }
        }

        /*
         * Raven
         */
        if (!params.skip_raven) {
          raven_assembly(lreads)
          pilon_polish_raven(raven_assembly.out[1], preads.concat(sreads).collect())
          raven_ch = pilon_polish_raven.out[1]

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_raven(raven_assembly.out[1], lreads, fast5, fast5_dir)
            pilon_polish_raven_nanopolish(nanopolish_raven.out[0], preads.concat(sreads).collect())
            nanopol_raven_ch = pilon_polish_raven_nanopolish.out[1]
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_raven(raven_assembly.out[1], lreads)
            pilon_polish_raven_medaka(medaka_raven.out[1], preads.concat(sreads).collect())
            medaka_raven_ch = pilon_polish_raven_medaka.out[1]
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_raven(raven_assembly.out[1], bamFile, nBams)
            pilon_polish_raven_variantCaller(variantCaller_raven.out[1], preads.concat(sreads).collect())
            arrow_raven_ch = pilon_polish_raven_variantCaller.out[1]
          }
        }
      }

      // Run quast
      quast(spades_ch.mix(spades_ch, unicycler_hybrid_ch, haslr_ch,
                          canu_ch, nanopol_canu_ch, medaka_canu_ch, arrow_canu_ch,
                          raven_ch, nanopol_raven_ch, medaka_raven_ch, arrow_raven_ch,
                          flye_ch, nanopol_flye_ch, medaka_flye_ch, arrow_flye_ch,
                          unicycler_ch, nanopol_unicycler_ch, medaka_unicycler_ch, arrow_unicycler_ch), preads.concat(sreads).collect())

      // Run multiqc
      multiqc(quast.out[1].collect(), Channel.value('hybrid'))

}
