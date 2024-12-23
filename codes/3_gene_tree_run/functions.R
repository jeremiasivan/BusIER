# functions for codes/3_gene_tree_run

# function: quality control for short reads
f_qc_short_reads <- function(fastq_one, fastq_two, fn_adapters, prefix, min_quality, thread, exe_adapterremoval) {
    cmd_qc <- paste(exe_adapterremoval, "--file1", fastq_one)
    
    # check the reverse fastq file
    if (!is.null(fastq_two)) {
        cmd_qc <- paste(cmd_qc, "--file2", fastq_two)
    }

    # update the adapters
    if (fn_adapters != "" && file.exists(fn_adapters)) {
        cmd_qc <- paste(cmd_qc, "--adapter-list", fn_adapters)
    }

    # set the minimum quality score and merge overlapping reads
    cmd_qc <- paste(cmd_qc, "--basename", prefix,
                    "--trimqualities --trimns --minquality", min_quality)
    
    # run AdapterRemoval to get prefix.collapsed.truncated
    system(cmd_qc)
}

# function: build database for Easy353
f_easy353_build_database <- function(taxonomy, dir_out, thread, exe_build_db) {
    cmd_build_db <- paste(exe_build_db,
                          "-o", dir_out,
                          "-c", taxonomy,
                          "-t", thread, "-generate")
    system(cmd_build_db)
}

# function: run Easy353
f_easy353_run <- function(fastq_one, fastq_two, dir_db, dir_out, thread, exe_easy353) {
    cmd_easy353 <- paste(exe_easy353, "-1", fastq_one)

    if (!is.null(fastq_two)) {
        cmd_easy353 <- paste(cmd_easy353, "-2", fastq_two)
    }

    cmd_easy353 <- paste(cmd_easy353, "-r", dir_db, "-o", dir_out, "-ft 1 -at", thread)
    system(cmd_easy353)
}