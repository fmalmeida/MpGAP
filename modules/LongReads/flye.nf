process flye_assembly {
  publishDir "${params.outdir}/${lrID}/longreads_only", mode: 'copy', overwrite: true
  container 'fmalmeida/mpgap'
  cpus params.threads
  tag "Performing a longreads only assembly with Flye"

  input:
  file lreads

  output:
  file "flye" // Saves all files
  tuple file("flye/assembly_flye.fasta"), val(lrID), val('flye') // Gets contigs file

  script:
  lr    = (params.lr_type == 'nanopore') ? '--nano-raw' : '--pacbio-raw'
  lrID  = lreads.getSimpleName()
  gsize = (params.genomeSize) ? "--genome-size ${params.genomeSize}" : ""
  """
  source activate flye ;
  flye ${lr} $lreads ${gsize} --out-dir flye \
  --threads ${params.threads} ${params.flye_additional_parameters} &> flye.log ;
  mv flye/assembly.fasta flye/assembly_flye.fasta
  """
}
