process longstich {
    publishDir "${params.output}/${prefix}/longstich_scaffolded", mode: 'copy'
    tag "${id}"
    label 'process_assembly'

    input:
    tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

    output:
    file "${assembler}_gcpp_variants.gff" // Save gff
    tuple val(id), file("${assembler}_gcpp_consensus.fasta"), val("${assembler}_gcpp") // Save contigs
    path('versions.yml'), emit: versions

    when:
    entrypoint == 'hybrid_strategy_3'

    script:
    def version   = '1.0.5' // does not have version command
    if ( lr_type == 'nanopore' ) { reads_map = 'ont' }
    if ( lr_type == 'pacbio'   ) { reads_map = high_quality_longreads ? 'hifi' : 'pb' }
    """
    # run longstich
    longstitch \\
        run \\
        draft=${draft} \\
        reads=${lreads} \\
        G=${genome_size} \\
        t=${task.cpus} \\
        out_prefix=${id} \\
        rounds=4 \\
        longmap=${reads_map}

    # get version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longstich: ${version}
    END_VERSIONS
    """
}
