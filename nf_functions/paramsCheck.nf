
def paramsCheck() {

  /*
   * Checking for nf tower parameters
   */
  if ((params.use_tower && !params.tower_token) || (!params.use_tower && params.tower_token)) {
    println """
    ERROR!
    A minor error has occurred
      ==> User forgot to set (together) both --use_tower and --tower_token. Or they are mispelled.
    These parameters must be used together. They are the necessary to run the pipeline using the amazing Nextflow Tower!
    Cheers.
    """.stripIndent()

    exit 1
  }

  /*
   * Checking for genomeSize option when using canu
   */
  if (!params.skip_canu && !params.genomeSize && !params.shortreads_paired && !params.shortreads_single
      && params.longreads && params.lr_type) {
    println """
    ERROR!
    A minor error has occurred
      ==> The pipeline will try to run canu assembler, but user forgot to tell the expected genome size.
    Please either set the expected genome size with --genomeSize, or turn off the canu assembly with --skip_canu
    Cheers.
    """.stripIndent()

    exit 1
  }

  /*
   * Checking for genomeSize option when using haslr
   */
  if (!params.skip_haslr && !params.genomeSize && (params.shortreads_paired || params.shortreads_single)
      && params.longreads && params.lr_type) {
    println """
    ERROR!
    A minor error has occurred
      ==> The pipeline will try to run haslr assembler, but user forgot to tell the expected genome size.
    Please either set the expected genome size with --genomeSize, or turn off the haslr assembly with --skip_haslr
    Cheers.
    """.stripIndent()

    exit 1
  }

}
