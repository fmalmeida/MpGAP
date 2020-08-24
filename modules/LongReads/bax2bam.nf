process bax2bam {
  publishDir "${params.outdir}/longreads-only/${params.prefix}_subreads", mode: 'copy', overwrite: true
  container 'fmalmeida/mpgap'
  cpus params.threads

  input:
  file(bax)

  output:
  file("*.subreads.bam") // Get all bam files produced

  script:
  """
  source activate pacbio ;
  bax2bam ${bax.join(" ")} --subread  \
  --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,MergeQV,SubstitutionQV,PulseWidth,SubstitutionTag;
  """
}
