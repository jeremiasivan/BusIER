# functions for codes/2_busco_check

# function: read map fastq to reference sequence using BWA-MEM2
f_read_mapping <- function(refseq, fastq, thread, exe_bwamem2, file_sam) {
    # index reference file
    cmd_index <- paste(exe_bwamem2, "index", refseq)
    system(cmd_index)

    # read-map
    cmd_readmap <- paste(exe_bwamem2, "mem -t", thread, refseq)

    len_fastq <- length(fastq)
    if (len_fastq == 1) {
        cmd_readmap <- paste(cmd_readmap, fastq[1], ">", file_sam)
    } else if (len_fastq == 2) {
        cmd_readmap <- paste(cmd_readmap, fastq[1], fastq[2], ">", file_sam)
    }
    system(cmd_readmap)
}

# function: convert SAM to BAM and VCF
f_variant_calling <- function(prefix, dir_output, thread, refseq, exe_samtools, exe_bcftools) {
    # initiate variables
    fn_sam <- paste0(dir_output, "/", prefix, ".sam")
    fn_bam <- paste0(dir_output, "/", prefix, ".bam")
    fn_vcf <- paste0(dir_output, "/", prefix, ".vcf.gz")

    nthread <- paste("--threads", thread)

    # run samtools fixmate
    cmd_samtools <- paste(exe_samtools, "view", nthread, "-b -u", fn_sam, "|",
                          exe_samtools, "collate", nthread, "-O -u - |",
                          exe_samtools, "fixmate", nthread, "-m -u - - |",
                          exe_samtools, "sort", nthread, "-u - |",
                          exe_samtools, "markdup", nthread, "-", fn_bam)
    system(cmd_samtools)                  
    
    # run bcftools mpileup
    cmd_bcftools <- paste(exe_bcftools, "mpileup", nthread, "-Ou -f", refseq, fn_bam, "|",
                          exe_bcftools, "call", nthread, "-Ou -mv |",
                          exe_bcftools, "norm", nthread, "-f", refseq, "-Oz -o", fn_vcf) 
    system(cmd_markdup)
}

# function: run BUSCO pipeline
f_run_busco <- function(fn_fasta, lineage, prefix, dir_output, mode, thread, exe_busco) {
    cmd_busco <- paste(exe_busco,
                       "--offline",
                       "-i", fn_fasta,
                       "-l", lineage,
                       "-m", mode,
                       "-o", prefix,
                       "--out_path", dir_output,
                       "-c", thread,
                       "--quiet --force")
    system(cmd_busco)
}

# functions: extract BUSCO region
f_extract_busco <- function(fn_vcf, refseq, coordinates, fn_out, exe_samtools, exe_bcftools) {
    # extract consensus BUSCO sequence
    cmd_consensus <- paste(exe_samtools, "faidx", refseq, coordinates, "|",
                           exe_bcftools, "consensus", fn_vcf, "-o", fn_out)                                    
    system(cmd_consensus)
}

# functions: combine individual FASTA as MSA
f_fasta2msa <- function(fn_input, header, fn_out) {
    # initiate variable
    first_sequence <- TRUE

    # open the FASTA file
    con <- file(fn_input, "r")

    # iterate over lines
    while (length(line <- readLines(con, n = 1)) > 0) {
        if (grepl("^>+", line)) {
            if (first_sequence) {
                write.table(paste0(">", header), file=fn_out, quote=F, row.names=F, col.names=F, append=T)
                first_sequence <- FALSE
            }
        } else {
            write.table(line, file=fn_out, quote=F, row.names=F, col.names=F, append=T)
        }
    }
    
    # close the file connection
    close(con)
}

# functions: run IQ-Tree 2
f_iqtree2 <- function(fn_input, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-T 1 --quiet -redo")
    system(cmd_iqtree2)
}