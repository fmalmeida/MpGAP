.. _examples:

CLI usage Examples
==================

Illumina paired end reads.
""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir outputs/illumina_paired --run_shortreads_pipeline --shortreads
      "illumina/SRR9847694_{1,2}.fastq.gz" --reads_size 2 --lighter_execute --lighter_genomeSize 4600000 --clip_r1 5 --three_prime_clip_r1 5
      --clip_r2 5 --three_prime_clip_r2 5 --quality_trim 30 --flash_execute

.. note::

  Since it will always be a pattern match, example "illumina/SRR9847694_{1,2}.fastq.gz", it MUST ALWAYS be double quoted as the example below.

Illumina single end reads.
""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/illumina_single --run_shortreads_pipeline
      --shortreads "sample_dataset/illumina/SRR9696*.fastq.gz" --reads_size 1 --clip_r1 5 --three_prime_clip_r1 5

.. note::

  Multiple files at once, using fixed number of bases to be trimmed. If multiple unpaired reads are given as input at once, pattern MUST be double quoted: "SRR9696*.fastq.gz"

ONT reads (fastq)
"""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/ont --run_longreads_pipeline --lreads_type nanopore
  --longReads sample_dataset/ont/kpneumoniae_25X.fastq --nanopore_prefix kpneumoniae_25X

Pacbio basecalled (.fastq) reads with nextflow general report
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/pacbio_from_fastq --run_longreads_pipeline --lreads_type pacbio
  --longReads sample_dataset/pacbio/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreads.fastq -with-report

Pacbio raw (subreads.bam) reads
"""""""""""""""""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/pacbio --run_longreads_pipeline --lreads_type pacbio
  --pacbio_bamPath sample_dataset/pacbio/m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.subreads.bam

Pacbio raw (legacy .bas.h5 to subreads.bam) reads
"""""""""""""""""""""""""""""""""""""""""""""""""

::

  ./nextflow run fmalmeida/ngs-preprocess --threads 3 --outDir sample_dataset/outputs/pacbio --run_longreads_pipeline --lreads_type pacbio
  --pacbio_h5Path sample_dataset/pacbio/m140912_020930_00114_c100702482550000001823141103261590_s1_p0.1.bas.h5


Running with a configuration file
"""""""""""""""""""""""""""""""""

::

      ./nextflow run fmalmeida/ngs-preprocess -c nextflow.config
