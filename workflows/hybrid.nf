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
include { quast as quast_lreads_canu;       quast as quast_lreads_flye; quast as quast_lreads_unicycler; quast as quast_lreads_raven;
          quast as quast_nanopolish_canu;   quast as quast_nanopolish_flye; quast as quast_nanopolish_unicycler; quast as quast_nanopolish_raven;
          quast as quast_medaka_canu;       quast as quast_medaka_flye; quast as quast_medaka_unicycler; quast as quast_medaka_raven;
          quast as quast_variantcaller_canu; quast as quast_variantcaller_flye; quast as quast_variantcaller_unicycler;
          quast as quast_hybrid_unicycler; quast as quast_hybrid_haslr; quast as quast_hybrid_spades } \
       from '../modules/QualityAssessment/quast.nf' params(threads: params.threads, outdir: params.outdir, longreads: params.longreads,
            shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)

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
       * Full (default) hybrid mode
       */

      if (!params.strategy_2) {
        // SPAdes
        if (!params.skip_spades) {
          spades_hybrid(lreads, preads, sreads)
          quast_hybrid_spades(spades_hybrid.out[1], preads.concat(sreads).collect())
        }
        // Unicycler
        if (!params.skip_unicycler) {
          unicycler_hybrid(lreads, preads, sreads)
          quast_hybrid_unicycler(unicycler_hybrid.out[1], preads.concat(sreads).collect())
        }
        // Haslr
        if (!params.skip_haslr) {
          haslr_hybrid(lreads, preads, sreads)
          quast_hybrid_haslr(haslr_hybrid.out[1], preads.concat(sreads).collect())
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
          quast_lreads_canu(pilon_polish_canu.out[1], preads.concat(sreads).collect())

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_canu(canu_assembly.out[1], lreads, fast5, fast5_dir)
            pilon_polish_canu_nanopolish(nanopolish_canu.out[0], preads.concat(sreads).collect())
            quast_nanopolish_canu(pilon_polish_canu_nanopolish.out[1], preads.concat(sreads).collect())
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_canu(canu_assembly.out[1], lreads)
            pilon_polish_canu_medaka(medaka_canu.out[1], preads.concat(sreads).collect())
            quast_medaka_canu(pilon_polish_canu_medaka.out[1], preads.concat(sreads).collect())
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_canu(canu_assembly.out[1], bamFile, nBams)
            pilon_polish_canu_variantCaller(variantCaller_canu.out[1], preads.concat(sreads).collect())
            quast_variantcaller_canu(pilon_polish_canu_variantCaller.out[1], preads.concat(sreads).collect())
          }
        }

        /*
         * Flye
         */
        if (!params.skip_flye) {
          flye_assembly(lreads)
          pilon_polish_flye(flye_assembly.out[1], preads.concat(sreads).collect())
          quast_lreads_flye(pilon_polish_flye.out[1], preads.concat(sreads).collect())

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_flye(flye_assembly.out[1], lreads, fast5, fast5_dir)
            pilon_polish_flye_nanopolish(nanopolish_flye.out[0], preads.concat(sreads).collect())
            quast_nanopolish_flye(pilon_polish_flye_nanopolish.out[1], preads.concat(sreads).collect())
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_flye(flye_assembly.out[1], lreads)
            pilon_polish_flye_medaka(medaka_flye.out[1], preads.concat(sreads).collect())
            quast_medaka_flye(pilon_polish_flye_medaka.out[1], preads.concat(sreads).collect())
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_flye(flye_assembly.out[1], bamFile, nBams)
            pilon_polish_flye_variantCaller(variantCaller_flye.out[1], preads.concat(sreads).collect())
            quast_variantcaller_flye(pilon_polish_flye_variantCaller.out[1], preads.concat(sreads).collect())
          }

        }

        /*
         * Unicycler
         */
        if (!params.skip_unicycler) {
          unicycler_lreads(lreads)
          pilon_polish_unicycler(unicycler_lreads.out[1], preads.concat(sreads).collect())
          quast_lreads_unicycler(pilon_polish_unicycler.out[1], preads.concat(sreads).collect())

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_unicycler(unicycler_lreads.out[1], lreads, fast5, fast5_dir)
            pilon_polish_unicycler_nanopolish(nanopolish_unicycler.out[0], preads.concat(sreads).collect())
            quast_nanopolish_unicycler(pilon_polish_unicycler_nanopolish.out[1], preads.concat(sreads).collect())
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_unicycler(unicycler_lreads.out[1], lreads)
            pilon_polish_unicycler_medaka(medaka_unicycler.out[1], preads.concat(sreads).collect())
            quast_medaka_unicycler(pilon_polish_unicycler_medaka.out[1], preads.concat(sreads).collect())
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_unicycler(unicycler_lreads.out[1], bamFile, nBams)
            pilon_polish_unicycler_variantCaller(variantCaller_unicycler.out[1], preads.concat(sreads).collect())
            quast_variantcaller_unicycler(pilon_polish_unicycler_variantCaller.out[1], preads.concat(sreads).collect())
          }
        }

        /*
         * Raven
         */
        if (!params.skip_raven) {
          raven_assembly(lreads)
          pilon_polish_raven(raven_assembly.out[1], preads.concat(sreads).collect())
          quast_lreads_raven(pilon_polish_raven.out[1], preads.concat(sreads).collect())

          // Nanopolish?
          if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
            nanopolish_raven(raven_assembly.out[1], lreads, fast5, fast5_dir)
            pilon_polish_raven_nanopolish(nanopolish_raven.out[0], preads.concat(sreads).collect())
            quast_nanopolish_raven(pilon_polish_raven_nanopolish.out[1], preads.concat(sreads).collect())
          }

          // Medaka?
          if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
            medaka_raven(raven_assembly.out[1], lreads)
            pilon_polish_raven_medaka(medaka_raven.out[1], preads.concat(sreads).collect())
            quast_medaka_raven(pilon_polish_raven_medaka.out[1], preads.concat(sreads).collect())
          }

          // VariantCaller?
          if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
            variantCaller_raven(raven_assembly.out[1], bamFile, nBams)
            pilon_polish_raven_variantCaller(variantCaller_raven.out[1], preads.concat(sreads).collect())
            quast_variantcaller_raven(pilon_polish_raven_variantCaller.out[1], preads.concat(sreads).collect())
          }
        }
      }

}
