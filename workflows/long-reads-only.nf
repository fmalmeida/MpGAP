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
 * Module for assessing assembly qualities
 */
include { quast as quast_lreads_canu;       quast as quast_lreads_flye; quast as quast_lreads_unicycler; quast as quast_lreads_raven;
          quast as quast_nanopolish_canu;   quast as quast_nanopolish_flye; quast as quast_nanopolish_unicycler; quast as quast_nanopolish_raven;
          quast as quast_medaka_canu;       quast as quast_medaka_flye; quast as quast_medaka_unicycler; quast as quast_medaka_raven;
          quast as quast_variantcaller_canu; quast as quast_variantcaller_flye; quast as quast_variantcaller_unicycler } \
       from '../modules/QualityAssessment/quast.nf' params(threads: params.threads, outdir: params.outdir, longreads: params.longreads,
            shortreads_paired: params.shortreads_paired, shortreads_single: params.shortreads_single, lr_type: params.lr_type)

workflow lreadsonly_nf {
  take:
      reads
      fast5
      fast5_dir
      bamFile
      nBams
  main:
      /*
       * Canu
       */
      if (!params.skip_canu) {
        canu_assembly(reads)
        quast_lreads_canu(canu_assembly.out[1], reads)

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_canu(canu_assembly.out[1], reads, fast5, fast5_dir)
          quast_nanopolish_canu(nanopolish_canu.out[0], reads)
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_canu(canu_assembly.out[1], reads)
          quast_medaka_canu(medaka_canu.out[1], reads)
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_canu(canu_assembly.out[1], bamFile, nBams)
          quast_variantcaller_canu(variantCaller_canu.out[1], reads)
        }
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        flye_assembly(reads)
        quast_lreads_flye(flye_assembly.out[1], reads)

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_flye(flye_assembly.out[1], reads, fast5, fast5_dir)
          quast_nanopolish_flye(nanopolish_flye.out[0], reads)
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_flye(flye_assembly.out[1], reads)
          quast_medaka_flye(medaka_flye.out[1], reads)
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_flye(flye_assembly.out[1], bamFile, nBams)
          quast_variantcaller_flye(variantCaller_flye.out[1], reads)
        }

      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        unicycler_lreads(reads)
        quast_lreads_unicycler(unicycler_lreads.out[1], reads)

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_unicycler(unicycler_lreads.out[1], reads, fast5, fast5_dir)
          quast_nanopolish_unicycler(nanopolish_unicycler.out[0], reads)
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_unicycler(unicycler_lreads.out[1], reads)
          quast_medaka_unicycler(medaka_unicycler.out[1], reads)
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_unicycler(unicycler_lreads.out[1], bamFile, nBams)
          quast_variantcaller_unicycler(variantCaller_unicycler.out[1], reads)
        }
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        raven_assembly(reads)
        quast_lreads_raven(raven_assembly.out[1], reads)

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_raven(raven_assembly.out[1], reads, fast5, fast5_dir)
          quast_nanopolish_raven(nanopolish_raven.out[0], reads)
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_raven(raven_assembly.out[1], reads)
          quast_medaka_raven(medaka_raven.out[1], reads)
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_raven(raven_assembly.out[1], bamFile, nBams)
          quast_variantcaller_raven(variantCaller_raven.out[1], reads)
        }
      }
}
