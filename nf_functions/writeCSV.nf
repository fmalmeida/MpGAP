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
      it.pacbio_bam = (it.pacbio_bam) ? it.pacbio_bam : "missing_pacbio_bam"

      /*
       * Check for wtdbg2 technology
       */
      it.wtdbg2_technology = (it.wtdbg2_technology) ? it.wtdbg2_technology : params.wtdbg2_technology

      /*
       * Create entrypoints
       */

       // has only short reads
       if ((it.paired != "missing_paired" || it.single != "missing_single") && it.lreads == "missing_lreads") {
         entrypoint = 'shortreads_only'
       } 
       
       // has both short reads and long reads
       else if (((it.paired != "missing_paired" || it.single != "missing_single") && it.lreads != "missing_lreads")) {
         if (it.hybrid_strategy) {
           entrypoint = 'hybrid_strategy_' + it.hybrid_strategy.toString().toLowerCase()
         } else if (params.strategy_2) {
           entrypoint = 'hybrid_strategy_2'
         } else {
           entrypoint = 'hybrid_strategy_1'
         }
         if (entrypoint != 'hybrid_strategy_1' && entrypoint != 'hybrid_strategy_2' && entrypoint != 'hybrid_strategy_both') {
           println "ERROR: In the YAML, the 'strategy:' key must be either 1, 2 or both"
           exit 1
         }
       }
       
       // has only long reads
       else if (((it.paired == "missing_paired" || it.single == "missing_single") && it.lreads != "missing_lreads")) {
         entrypoint = 'longreads_only'
       } 
       
       // has nothing
       else {
         println "ERROR: At least one read type must be set: illumina, nanopore or pacbio!"
         exit 1
       }
      
      // output csv
      if (entrypoint == 'hybrid_strategy_both') {
        // hybrid strategy 1 and 2
        // separated with \n
        "${it.id}:strategy_1,hybrid_strategy_1,${it.fwd_pair},${it.rev_pair},${it.single},${it.lreads},${it.lr_type},${it.wtdbg2_technology},${it.genomeSize},${it.corrected_lreads},${it.medaka_model},${it.fast5},${it.pacbio_bam}\n${it.id}:strategy_2,hybrid_strategy_2,${it.fwd_pair},${it.rev_pair},${it.single},${it.lreads},${it.lr_type},${it.wtdbg2_technology},${it.genomeSize},${it.corrected_lreads},${it.medaka_model},${it.fast5},${it.pacbio_bam}"
      } else {
        "${it.id},${entrypoint},${it.fwd_pair},${it.rev_pair},${it.single},${it.lreads},${it.lr_type},${it.wtdbg2_technology},${it.genomeSize},${it.corrected_lreads},${it.medaka_model},${it.fast5},${it.pacbio_bam}"
      }
  }
}
