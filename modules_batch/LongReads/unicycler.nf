process unicycler {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "${id}: unicycler assembly"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model),file(fast5), val(fast5_dir), file(bams), val(nBams), val(prefix)

  output:
  file "unicycler/*" // Save all files
  tuple val(id), file("unicycler/unicycler_assembly.fasta"), val('unicycler') // Gets contigs file

  script:
  """
  unicycler -l $lreads \
  -o unicycler -t ${params.threads} \
  ${params.unicycler_additional_parameters}

  # Rename assembly
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}
