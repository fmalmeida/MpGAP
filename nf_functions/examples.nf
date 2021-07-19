def exampleMessage() {
    log.info """

    Examplification on how to run fmalmeida/mpgap pipeline using the CLI configuration

    ## Short reads only - PAIRED:

\$ nextflow run fmalmeida/mpgap --outdir output --threads 5 --shortreads_paired "path-to/illumina_r{1,2}.fastq"

    ## Short reads only - SINGLE:

\$ nextflow run fmalmeida/mpgap --outdir output --threads 5 --shortreads_single "path-to/illumina_unpaired.fastq"

    ## Short reads only - Both PAIRED and SINGLE:

\$ nextflow run fmalmeida/mpgap --outdir output --threads 5 --shortreads_paired "path-to/illumina_r{1,2}.fastq" --shortreads_single "path-to/illumina_unpaired.fastq"

    ## Long reads only - ONT (Using both Nanopolish and Medaka):

\$ nextflow run fmalmeida/mpgap --outdir output --threads 5 --medaka_sequencing_model r941_min_fast_g303 --nanopolish_fast5Path "path-to/fast5_pass" \
--nanopolish_max_haplotypes 2000 --genomeSize 2m --lr_type nanopore --longreads "path-to/ont_reads.fastq"

    ## Long reads only - Pacbio (Using gcpp):

\$ nextflow run fmalmeida/mpgap --outdir output --threads 5 --genomeSize 2m --lr_type "pacbio" --longreads "path-to/pacbio.subreads.fastq" --pacbio_bams "path-to/pacbio.*.subreads.bam"

    ## Hybrid assembly - Using both paired and single end short reads:

\$ nextflow run fmalmeida/mpgap --outdir output --threads 5 --shortreads_paired "path-to/illumina_r{1,2}.fastq" \
--shortreads_single "path-to/illumina_unpaired.fastq" --lr_type nanopore --longreads "path-to/ont_reads.fastq" --genomeSize 4m

    ## Hybrid assembly - by polishing (with shortreads) a longreads-only assembly
    ## Also using medaka prior to the polishing with shortreads

\$ nextflow run fmalmeida/mpgap --outdir output --threads 5 --shortreads_paired "path-to/illumina_r{1,2}.fastq" \
--lr_type nanopore --longreads "path-to/ont_reads.fastq" --strategy_2 --medaka_sequencing_model r941_min_fast_g303 --genomeSize 4m
    """.stripIndent()
}
