def exampleMessage() {
    log.info """

    Examplification on how to run fmalmeida/MpGAP pipeline using the CLI configuration

    Short reads only - PAIRED:
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type illumina-only --try_spades --try_unicycler --shortreads_paired "dataset_1/sampled/illumina_R{1,2}.fastq"

    Short reads only - SINGLE:
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type illumina-only --try_spades --try_unicycler --shortreads_single "dataset_1/sampled/illumina_single.fastq"

    Short reads only - Both PAIRED and SINGLE:
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type illumina-only --try_spades --try_unicycler --shortreads_paired "dataset_1/sampled/illumina_R{1,2}.fastq" --shortreads_single "dataset_1/sampled/illumina_single.fastq"

    Long reads only - ONT (Using both Nanopolish and Medaka):
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type longreads-only --try_canu --try_flye --try_unicycler --medaka_sequencing_model r941_min_fast_g303 \
--nanopolish_fast5Path "dataset_1/ont/fast5_pass" --nanopolish_max_haplotypes 2000 --genomeSize 2m --lr_type nanopore --longreads "dataset_1/ont/ont_reads.fastq"

    Long reads only - Pacbio (Using VariantCaller):
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type "longreads-only" --try_unicycler --try_flye --genomeSize "2m" --lr_type "pacbio" \
--longreads "E01_1/Analysis_Results/preprocessed/longreads/pacbio/m141013_011508_sherri_c100709962550000001823135904221533_s1_p0.subreads.subset.fastq" --pacbio_all_bam_path "E01_1/Analysis_Results/preprocessed/longreads/pacbio/*.subreads.bam"

    Hybrid assembly - Using both paired and single end short reads
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type hybrid --try_unicycler --shortreads_paired "dataset_1/sampled/illumina_R{1,2}.fastq" \
--shortreads_single "dataset_1/sampled/illumina_single.fastq" --lr_type nanopore --longreads "dataset_1/ont/ont_reads.fastq"

    Hybrid assembly - by polishing (with shortreads) a longreads-only assembly -- also using medaka
\$ nextflow run fmalmeida/MpGAP --outdir output --threads 5 --assembly_type hybrid --try_unicycler --try_flye --try_canu --shortreads_paired "dataset_1/sampled/illumina_R{1,2}.fastq"
--genomeSize 2m --lr_type nanopore --longreads "dataset_1/ont/ont_reads.fastq" --illumina_polish_longreads_contigs --medaka_sequencing_model r941_min_fast_g303
    """.stripIndent()
}
