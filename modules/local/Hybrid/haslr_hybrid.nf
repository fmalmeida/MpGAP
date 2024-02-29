process haslr_hybrid {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "*" // Save everything
  tuple val(id), file("haslr/haslr_assembly.fa"), val('haslr') // Gets contigs file
  path('versions.yml'), emit: versions

  when:
  ((!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) || !(single =~ /input.*/)) && !(lreads =~ /input.*/) && (entrypoint == 'hybrid_strategy_1')

  script:
  // Check reads
  paired_reads = (!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) ? "$sread1 $sread2" : ""
  single_reads = !(single =~ /input.*/) ? "$single" : ""
  additional_params = (params.haslr_additional_parameters) ? params.haslr_additional_parameters : ""
  """
  # run haslr
  /opt/haslr/bin/haslr.py \\
      -t $task.cpus \\
      -o haslr \\
      -g ${genome_size} \\
      -l $lreads \\
      -x ${lr_type} \\
      $additional_params \\
      -s ${paired_reads} ${single_reads} 

  # rename results
  if [[ -f haslr/*/asm.final.fa ]]
  then
    cp haslr/*/asm.final.fa haslr/haslr_assembly.fa
  else
    echo 'Haslr did not generate any assembly. This may be due the reads, but also can be due haslr tool itself.'
    echo 'See this issue: https://github.com/vpc-ccg/haslr/issues/4'
    echo 'If your you know your reads are good, then consider running the pipeline without haslr, --skip_haslr.'
    exit 2
  fi

  # get version
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      haslr: \$( haslr.py --version )
  END_VERSIONS
  """
}
