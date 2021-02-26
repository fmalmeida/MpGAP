
def paramsCheck() {

  /*
   * Checking for genomeSize option when using canu
   */
  if (!params.skip_canu && !params.genomeSize && params.longreads &&
    params.lr_type && !params.shortreads_paired && !params.shortreads_single) {
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
   * Checking for genomeSize option when using canu
   */
  if (!params.skip_canu && !params.genomeSize && params.longreads && params.strategy_2 &&
    params.lr_type && (params.shortreads_paired || params.shortreads_single)) {
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
  if (!params.skip_haslr && !params.genomeSize
      && (params.shortreads_paired || params.shortreads_single)
      && params.longreads && params.lr_type && !params.strategy_2) {
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
