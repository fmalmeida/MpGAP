process medaka {
    publishDir "${params.output}/${prefix}/medaka_polished_contigs", mode: 'copy'
    tag "${id}"
    label 'process_assembly'

    input:
    tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

    output:
    file "${assembler}" // Save everything
    tuple val(id), file("${assembler}/${assembler}_medaka_consensus.fa"), val("${assembler}_medaka") // Save medaka contigs
    path('versions.yml'), emit: versions

    when:
    (medaka_model) && (lr_type == 'nanopore') && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

    script:
    """
    # map reads
    minimap2 \
        -x map-ont \\
        ${draft} \\
        ${lreads} > reads_mapped.paf ;

    # first step racon polish
    # as in medaka manual
    racon \\
        -m 8 -x -6 -g -8 -w 500 \\
        -t $task.cpus \\
        ${lreads} \\
        reads_mapped.paf \\
        ${draft} > racon_consensus.fasta ;

    # second step medaka polish
    medaka_consensus \\
        -i ${lreads} \\
        -d racon_consensus.fasta \\
        -o ${assembler} \\
        -t $task.cpus \\
        -m ${medaka_model} ;

    # rename results
    mv ${assembler}/consensus.fasta ${assembler}/${assembler}_medaka_consensus.fa

    # get version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}
