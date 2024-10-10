# WDL Version Specification
version 1.0

# Task: mosdepthGCopt
# This WDL task runs the `mosdepth` tool to calculate the coverage over given intervals in CRAM files
# and performs GC normalization using an R script.

task mosdepthGCopt {
  input {
    Int threads = 4                       # Number of threads for `mosdepth` tool
    Int maxRetries = 0                     # Maximum number of retries for the task
    String batch_id = "1"                  # Identifier for the batch
    Array[File]+ cram_files                # List of input CRAM files
    Array[File]+ crai_files                # Corresponding CRAI index files
    File genome_fa                         # Reference genome in FASTA format
    File genome_fai                        # FASTA index file
    File genome_dict                       # Reference genome dictionary file
    File exome_capture_kit                 # BED file defining the exome capture regions
    File exome_capture_groups              # File defining exon groups for GC normalization
    Boolean fast_mode = false              # Fast mode flag (currently not used)
    File norm_Rscript                      # R script for GC content normalization
  }

  # Define the output file name based on batch_id
  String output_batch = batch_id + "_median_coverage"

  command <<<
    # Ensure strict error handling
    set -euo pipefail
    
    # Create output directory
    mkdir out_files
    
    # Loop through each CRAM file to calculate coverage
    for cram in ~{sep=" " cram_files}; do
      # Extract sample ID from the CRAM filename
      id=$(echo $cram | cut -d"/" -f6 | cut -d"." -f1 | cut -d"_" -f1)
      
      # Run `mosdepth` to calculate coverage using the exome capture kit
      mosdepth -n -t ~{threads} --by ~{exome_capture_kit} -m --fasta ~{genome_fa} $id $cram
      
      # Define the input file for the R script (output from `mosdepth`)
      input_R_file="${id}.regions.bed.gz"
      
      # Perform GC normalization using the R script
      Rscript ~{norm_Rscript} ${input_R_file} ~{exome_capture_groups}
    done
    
    # Change to the output directory
    cd out_files
    
    # Create a header for the final output file and merge all individual results
    echo "eid    chrY_norm" > ~{output_batch}
    cat [0-9]*.txt >> ~{output_batch}
    
    # Move the final output file to the main directory
    mv ~{output_batch} ../
  >>>

  # Define runtime settings
  runtime {
    docker: "./docker_images/chrycovergcnorm.tar.gz"  # Docker image for running the task
    maxRetries: maxRetries                                          # Set the maximum retries for the task
  }

  # Define the output file produced by the task
  output {
    File median_coverage = "~{output_batch}"                        # Output file with median coverage values
  }
}
