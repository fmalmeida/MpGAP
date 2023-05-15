process pilon_polish {

    publishDir "${params.output}/${prefix}/pilon_polished_contigs", mode: 'copy'
    tag "${id}"
    label 'process_assembly'

    input:
    tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

    output:
    file("${assembler}/*") // Get everything
    tuple val(id), file("${assembler}/${assembler}_pilon_consensus.fasta"), val("${assembler}_pilon_polished")

    script:
    paired_cmd = (sread1 =~ /input.*/) ? "" : "${sread1} ${sread2}"
    single_cmd = (single =~ /input.*/) ? "" : "${single}"
    fixed_id   = id - ":strategy_2"
    """
    # get tools path
    pilonjar_path=\$(find \$CONDA_PREFIX -name "pilon.jar" | grep "mpgap" | head -n 1)

    # Create the results dir
    rm -rf ${assembler} && mkdir -p ${assembler}; # make sure its fresh

    # starts with received draft
    export draft_to_polish=${draft}

    # how many rounds?
    for round in {1..${params.pilon_polish_rounds}} ; do

        echo "Running on \${draft_to_polish}!"

        # Index and align reads with bwa
        bwa index \${draft_to_polish} ;
        bwa mem \
            -M \
            -t $task.cpus \
            \${draft_to_polish} \
            ${paired_cmd} ${single_cmd} > ${fixed_id}_${assembler}_\${round}_aln.sam ;
        
        # convert alignment and index
        samtools \
            view \
            -bS \
            ${fixed_id}_${assembler}_\${round}_aln.sam | \
            samtools sort > ${fixed_id}_${assembler}_\${round}_aln.bam ;
        samtools index ${fixed_id}_${assembler}_\${round}_aln.bam ;

        # Execute pilon
        java \\
            -Xmx${params.pilon_memory_limit}G \\
            -jar \$pilonjar_path \\
            --genome \${draft_to_polish} \\
            --bam ${fixed_id}_${assembler}_\${round}_aln.bam \\
            --output ${assembler}_pilon_consensus_\${round} \\
            --outdir pilon_consensus_\${round} &> pilon_round_\${round}.log
        mv ${fixed_id}_${assembler}_\${round}_aln.bam pilon_round_\${round}.log pilon_consensus_\${round}

        # reset variable so next round has a different starting point
        export draft_to_polish=pilon_consensus_\${round}/${assembler}_pilon_consensus_\${round}.fasta ;
    
    done

    # save results
    ln -rs pilon_consensus_\${round} ${assembler}
    """

}
