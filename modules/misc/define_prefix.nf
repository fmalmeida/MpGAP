process define_prefix {

  tag "parsing inputs and generating prefix"

  input:
    file(lreads)
    tuple val(pair_id), file(pair_1), file(pair_2)
    file(sreads)

  output:
    env PREFIX

  script:
  
  /*
   * longreads only workflows
   */
  if (params.longreads && !params.shortreads_single && !params.shortreads_paired) {

    type    = "longreads_only"
    lrID    = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))
    out_dir = (params.prefix) ? "${params.prefix}/${type}" : "${lrID}/${type}"

  }

  /*
   * shortreads only workflows
   */
  
  // only single end reads
  if (!params.longreads && params.shortreads_single && !params.shortreads_paired) {

    type    = "shortreads_only"
    srID    = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    out_dir = (params.prefix) ? "${params.prefix}/${type}" : "${srID}/${type}"

  }

  // only paired end reads
  if (!params.longreads && !params.shortreads_single && params.shortreads_paired) {

    type    = "shortreads_only"
    prID    = pair_id
    out_dir = (params.prefix) ? "${params.prefix}/${type}" : "${prID}/${type}"

  }

  // both single and paired end reads
  if (!params.longreads && params.shortreads_single && params.shortreads_paired) {

    type    = "shortreads_only"
    srID    = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    prID    = pair_id
    out_dir = (params.prefix) ? "${params.prefix}/${type}" : "${prID}_and_${srID}/${type}"

  }

  /*
   * hybrid workflows
   */
  if (params.longreads && (params.shortreads_single || params.shortreads_paired)) {

    // check strategy
    type = (params.strategy_2) ? 'hybrid_strategy_2' : 'hybrid_strategy_1'

    // required long reads
    lrID    = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))

    // only single end reads
    if (params.shortreads_single && !params.shortreads_paired) {
      srID    = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
      out_dir = (params.prefix) ? "${params.prefix}/${type}" : "${lrID}_and_${srID}/${type}"
    }

    // only paired end reads
    if (!params.shortreads_single && params.shortreads_paired) {
      prID    = pair_id
      out_dir = (params.prefix) ? "${params.prefix}/${type}" : "${lrID}_and_${prID}/${type}"
    }

    // both single and paired end reads
    if (params.shortreads_single && params.shortreads_paired) {
      srID    = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
      prID    = pair_id
      out_dir = (params.prefix) ? "${params.prefix}/${type}" : "${lrID}_and_${prID}_and_${srID}/${type}"
    }
  }

  """
  PREFIX=${out_dir}
  """
}
