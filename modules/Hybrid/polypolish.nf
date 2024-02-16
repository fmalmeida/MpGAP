process polypolish {

    publishDir "${params.output}/${prefix}/polypolish_polished_contigs", mode: 'copy'
    tag "${id}"
    label 'process_assembly'

    input:
    tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

    output:
    file("${assembler}/*") // Get everything
    tuple val(id), file("${assembler}/${assembler}_polypolish_consensus.fasta"), val("${assembler}_polypolish_consensus")

    script:
    paired   = (sread2 =~ /input.*/) ? "false" : "true"
    fixed_id = id - ":strategy_2"
    """

    # Create the results dir
    rm -rf ${assembler} && mkdir -p ${assembler}; # make sure its fresh

    # Index and align reads with bwa
    bwa index ${draft} ;
    bwa mem \
        -M \
        -t $task.cpus \
        ${draft} \
        ${sread1} > ${fixed_id}_${assembler}_1_aln.sam ;
    if [ "${paired}" == "true" ] ; then
        bwa mem \
            -M \
            -t $task.cpus \
            ${draft} \
            ${sread2} > ${fixed_id}_${assembler}_2_aln.sam ;
    fi

    # Execute polypolish
    if [ "${paired}" == "true" ] ; then
        polypolish_insert_filter.py \
            --in1 ${fixed_id}_${assembler}_1_aln.sam \
            --in2 ${fixed_id}_${assembler}_2_aln.sam \
            --out1 filtered_1.sam \
            --out2 filtered_2.sam ;
    else
        polypolish_insert_filter.py \
            --in1 ${fixed_id}_${assembler}_1_aln.sam \
            --out1 filtered_1.sam ;
    fi
    polypolish ${draft} filtered*.sam > ${assembler}/${assembler}_polypolish_consensus.fasta
    """

}
