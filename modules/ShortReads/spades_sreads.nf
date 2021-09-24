process spades {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  tag { x }
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2), file(sreads), val(prefix)

  output:
  file "spades" // Save all output
  tuple file("spades/spades_assembly.fasta"), val(out_ids), val('spades') // Gets contigs file

  script:
  // Check reads
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads"
    x = "Performing a illumina-only assembly with SPAdes, using paired and single end reads"
    srId = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    out_ids = "${id}_and_${srId}"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads"
    id = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    x = "Performing a illumina-only assembly with SPAdes, using single end reads"
    out_ids = "${id}"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-1 $sread1 -2 $sread2"
    x = "Performing a illumina-only assembly with SPAdes, using paired end reads"
    out_ids = "${id}"
  }

  """
  spades.py -o spades -t ${params.threads} ${params.spades_additional_parameters} $parameter

  # Rename assembly
  mv spades/contigs.fasta spades/spades_assembly.fasta
  """
}

// batch mode
process spades_batch {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  tag "${id}: spades assembly"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), val(prefix)

  output:
  file "${assembler}" // Save all output
  tuple file("${assembler}/shovill_${assembler}_final.fasta"), val(id), val("shovill_${assembler}"), val(prefix), file(sread1), file(sread2), file(single) 

  when:
  (sread1 !=~ /input.?/ || sread2 !=~ /input.?/) && (single =~ /input.?/)

  script:
  """
  # Activate env
  source activate shovill ;

  # Run
  shovill \
    --outdir ${assembler} \
    --assembler ${assembler} \
    --R1 $sread1 --R2 $sread2 \
    --cpus ${params.threads} \
    --trim \
    ${params.shovill_additional_parameters}

  # Rename assembly
  mv ${assembler}/contigs.fa ${assembler}/shovill_${assembler}_final.fasta
  """
}
