def write_csv(in_list) {

  return in_list.collectFile( name:"samples.txt", sort: { it[0] }, newLine: true ) {

      /*
       * Check for illumina reads
       */
      if (it.illumina) {
        if (it.illumina.size() == 1) {
          it.paired   = "missing_paired"
          it.fwd_pair = "missing_pairFWD"
          it.rev_pair = "missing_pairREV"
          it.single   = it.illumina[0]
        } else if (it.illumina.size() == 2) {
          it.fwd_pair = it.illumina[0]
          it.rev_pair = it.illumina[1]
          it.single   = "missing_single"
        } else if (it.illumina.size() == 3) {
          it.fwd_pair = it.illumina[0]
          it.rev_pair = it.illumina[1]
          it.single   = it.illumina[2]
        }
      } else {
        it.paired   = "missing_paired"
        it.fwd_pair = "missing_pairFWD"
        it.rev_pair = "missing_pairREV"
        it.single   = "missing_single"
      }

      /*
       * Check for long reads
       */
      if (it.nanopore) {
        it.lr_type = "nanopore"
        it.lreads  = it.nanopore
      } else if (it.pacbio) {
        it.lr_type = "pacbio"
        it.lreads  = it.pacbio
      } else {
        it.lr_type = "missing_lr_type"
        it.lreads  = "missing_lreads"
      }
     
      /*
       * Check if long reads are corrected
       */
      it.corrected_lreads = (params.corrected_lreads || it.corrected_lreads) ? 'true' : 'false'

      /*
       * Check for fast5 directory
       */
      it.fast5 = (it.fast5) ? it.fast5 : "missing_fast5"

      /*
       * Check for medaka model
       */
      it.medaka_model = (it.medaka_model) ? it.medaka_model : params.medaka_sequencing_model

      /*
       * Check for sample genomeSize
       */
      if (it.genomeSize) {
        it.genomeSize = it.genomeSize
      } else if (params.genomeSize) {
        it.genomeSize = params.genomeSize
      } else {
        it.genomeSize = 'missing_genomeSize'
      }

      /*
       * Check for pacbio bams
       */
      it.pacbio_bams = (it.pacbio_bams) ? it.pacbio_bams : "missing_pacbio_bams"

      /*
       * Check for wtdbg2 technology
       */
      it.wtdbg2_technology = (it.wtdbg2_technology) ? it.wtdbg2_technology : params.wtdbg2_technology

      /*
       * Create entrypoints
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
      
      // output csv
      "${it.id},${entrypoint},${it.fwd_pair},${it.rev_pair},${it.single},${it.lreads},${it.lr_type},${it.wtdbg2_technology},${it.genomeSize},${it.corrected_lreads},${it.medaka_model},${it.fast5},${it.pacbio_bams}"
  }
}
