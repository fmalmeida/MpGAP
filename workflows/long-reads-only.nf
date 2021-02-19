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
       * Quast channels
       */
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
       * Canu
       */
      if (!params.skip_canu) {
        canu_assembly(reads)
        canu_ch = canu_assembly.out[1]

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_canu(canu_assembly.out[1], reads, fast5, fast5_dir)
          nanopol_canu_ch = nanopolish_canu.out[0]
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_canu(canu_assembly.out[1], reads)
          medaka_canu_ch = medaka_canu.out[1]
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_canu(canu_assembly.out[1], bamFile, nBams)
          arrow_canu_ch = variantCaller_canu.out[1]
        }
      }

      /*
       * Flye
       */
      if (!params.skip_flye) {
        flye_assembly(reads)
        flye_ch = flye_assembly.out[1]

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_flye(flye_assembly.out[1], reads, fast5, fast5_dir)
          nanopol_flye_ch = nanopolish_flye.out[0]
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_flye(flye_assembly.out[1], reads)
          medaka_flye_ch = medaka_flye.out[1]
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_flye(flye_assembly.out[1], bamFile, nBams)
          arrow_flye_ch = variantCaller_flye.out[1]
        }

      }

      /*
       * Unicycler
       */
      if (!params.skip_unicycler) {
        unicycler_lreads(reads)
        unicycler_ch = unicycler_lreads.out[1]

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_unicycler(unicycler_lreads.out[1], reads, fast5, fast5_dir)
          nanopol_unicycler_ch = nanopolish_unicycler.out[0]
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_unicycler(unicycler_lreads.out[1], reads)
          medaka_unicycler_ch = medaka_unicycler.out[1]
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_unicycler(unicycler_lreads.out[1], bamFile, nBams)
          arrow_unicycler_ch = variantCaller_unicycler.out[1]
        }
      }

      /*
       * Raven
       */
      if (!params.skip_raven) {
        raven_assembly(reads)
        raven_ch = raven_assembly.out[1]

        // Nanopolish?
        if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') {
          nanopolish_raven(raven_assembly.out[1], reads, fast5, fast5_dir)
          nanopol_canu_ch = nanopolish_raven.out[0]
        }

        // Medaka?
        if (params.medaka_sequencing_model && params.lr_type == 'nanopore') {
          medaka_raven(raven_assembly.out[1], reads)
          medaka_raven_ch = medaka_raven.out[1]
        }

        // VariantCaller?
        if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') {
          variantCaller_raven(raven_assembly.out[1], bamFile, nBams)
          arrow_raven_ch = variantCaller_raven.out[1]
        }
      }

      // Run quast
      quast(canu_ch.mix(nanopol_canu_ch, medaka_canu_ch, arrow_canu_ch,
                          raven_ch, nanopol_raven_ch, medaka_raven_ch, arrow_raven_ch,
                          flye_ch, nanopol_flye_ch, medaka_flye_ch, arrow_flye_ch,
                          unicycler_ch, nanopol_unicycler_ch, medaka_unicycler_ch, arrow_unicycler_ch), reads)

      // Run multiqc
      multiqc(quast.out[1].collect(), quast.out[2].distinct(), Channel.value('longreads-only'))
}
