process hifiasm {
  publishDir "${params.output}/${prefix}/hifiasm", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "hifiasm*" // Saves all files
  tuple val(id), file("hifiasm_assembly.fasta"), val('hifiasm') // Gets contigs file

  when:
  ( entrypoint == 'longreads_only'  || entrypoint == 'hybrid_strategy_2'  ) &&
  ( corrected_longreads.toBoolean() || high_quality_longreads.toBoolean() ) // only run hifiasm if reads are at least of good quality

  script:
  def additional_params = (params.hifiasm_additional_parameters) ? params.hifiasm_additional_parameters : ""
  """
  # run hifiasm
  hifiasm \\
      -o hifiasm \\
      -t $task.cpus \\
      $additional_params \\
      $lreads

  # convert to fasta
  awk '/^S/{print ">"\$2"\\n"\$3}' hifiasm.bp.p_ctg.gfa > hifiasm_assembly.fasta

  # check if assembly was produced
  if [ -s hifiasm_assembly.fasta ]
  then
    echo 'assembly: ok'
  else
    echo
    echo
    echo 'ERROR:'
    echo 'Assembly is empty. Are your reads really hifi?'
    echo 'Consider reviewing your reads or using --skip_hifiasm.'
    echo
    echo
    exit 2
  fi
  """
}
