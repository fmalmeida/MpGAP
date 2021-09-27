def write_csv(in_list) {

  return in_list.collectFile( name:"samples.txt", sort: { it[0] }, newLine: true ) {

    // Check short read pairs
    fwd_pair = (it.paired != "missing_paired") ? "${it.paired[0]}" : "missing_pairFWD"
    rev_pair = (it.paired != "missing_paired") ? "${it.paired[1]}" : "missing_pairREV"

    /*
     * Create entrypoints
     */

    /*
     * Short reads only
     */
    if ((it.paired != "missing_paired" || it.single != "missing_single") && it.lreads == "missing_lreads") {
      entrypoint = 'sr-only'
    } else if (((it.paired != "missing_paired" || it.single != "missing_single") && it.lreads != "missing_lreads")) {
      if (it.strategy == '2') {
        entrypoint = 'hybrid-strategy-2'
      } else {
        entrypoint = 'hybrid-strategy-1'
      }
    } else if (((it.paired == "missing_paired" || it.single == "missing_single") && it.lreads != "missing_lreads")) {
      entrypoint = 'lr-only'
    } else {
      println "ERROR: At least one read type must be set: illumina, nanopore or pacbio!"
      exit 1
    }

  /*
   * Output tuple
   */
  "${it.id},${entrypoint},${fwd_pair},${rev_pair},${it.single},${it.lreads},${it.lr_type},${it.fast5},${it.pacbio_bams}"
  }
}
