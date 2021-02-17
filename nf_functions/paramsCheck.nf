
def paramsCheck() {

  /*
      Checking for nf tower parameters
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

}
