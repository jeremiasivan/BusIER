# This script aims to extract intronic regions from BUSCOs. There should be only one reference sequence, but multiple mapped reads.
# Written by: Jeremias Ivan
# Last updated: 18 June 2024

###################################################################
######                     INPUT VARIABLES                   ######
###################################################################
file_metadata <- "/home/jeremias/BusIER/addons/busco_intron_metadata.tsv"

dir_busco_ref <- "/data/jeremias/run_eudicots_odb10/busco_sequences/single_copy_busco_sequences/"
dir_out <- "/data/jeremias/eucs/intron/"

exe_gffread <- "/home/jeremias/gffread/gffread"
exe_iqtree2 <- "/home/jeremias/iqtree-2.2.2.7-Linux/bin/iqtree2"
exe_astral <- "/home/jeremias/Astral/astral.5.7.8.jar"

thread <- 50

###################################################################
######                        FUNCTIONS                      ######
###################################################################

# function: add entry to logfile
f_write_log <- function(fn_log, msg) {
    write.table(msg, file=fn_log, quote=F, row.names=F, col.names=F, append=T)
}

# function: extract BUSCO coordinates from FASTA header
f_extract_coordinates <- function(fasta_header, busco, prefix) {
    # remove > sign
    no_header <- unlist(strsplit(fasta_header, split=">"))[2]
    ls_header <- unlist(strsplit(no_header, split=":"))
    ls_coordinates <- unlist(strsplit(ls_header[2], split="-"))

    seq_name <- ls_header[1]
    first_coordinate <- as.numeric(ls_coordinates[1]) + 1
    second_coordinate <- as.numeric(ls_coordinates[2]) + 1

    # set the start and stop coordinates
    strand <- "+"
    start_coordinate <- NULL
    stop_coordinate <- NULL

    if (first_coordinate < second_coordinate) {
        start_coordinate <- first_coordinate
        stop_coordinate <- second_coordinate
    } else if (first_coordinate > second_coordinate) {
        start_coordinate <- second_coordinate
        stop_coordinate <- first_coordinate
        strand <- "-"
    } else {
        return(list(errmsg=paste0("Error: ", busco, " coordinates for ", prefix, ". Skipped.")))
    }

    return(list(seqname=seq_name, start=start_coordinate, stop=stop_coordinate, strand=strand))
}

# function: manipulate and save GFF file
f_manipulate_gff <- function(fn_input, coordinates, busco, prefix, fn_out) {
    # check if entries comes from multiple genes
    f_check_target_id <- function(df_gff) {
        # extract attributes
        ls_attribute <- df_gff$attribute
        ls_target_id <- sapply(ls_attribute, function(x){
            stringr::str_match(x, "Target_ID=([^;]+)")[2]
        })

        if (length(unique(ls_target_id)) != 1) {
            # extract the most common target ID
            tb_target_id <- table(ls_target_id)
            target_id <- names(tb_target_id[tb_target_id==max(tb_target_id)])[1]

            # extract the GFF entries
            df_gff_filter <- df_gff[grepl(target_id, df_gff$attribute),]
            return(df_gff_filter)
        }

        return(df_gff)
    }

    # add intron entry to the gff
    f_extract_intron_gff <- function(first_exon, second_exon, i, busco, strand) {
        # extract coordinates
        start_coordinate <- first_exon$end + 1
        stop_coordinate <- second_exon$start - 1

        if (strand == "-") {
            start_coordinate <- second_exon$end + 1
            stop_coordinate <- first_exon$start - 1
        }

        return(data.table::data.table(seqname=first_exon$seqname,
                                      source=first_exon$source,
                                      feature="intron",
                                      start=start_coordinate,
                                      end=stop_coordinate,
                                      score=".",
                                      strand=strand,
                                      frame=".",
                                      attribute=paste0("ID=",busco,"_intron_",i)))
    }

    # read GFF file
    df_gff <- data.table::fread(fn_input, header=FALSE, select=1:9)
    data.table::setnames(df_gff, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

    # convert attribute to match with gffread namingconvention
    df_gff$attribute <- gsub("TCS_ID", "ID", df_gff$attribute)

    # extract relevant GFF entry
    df_gff <- df_gff[df_gff$seqname==coordinates$seqname & df_gff$strand==coordinates$strand]
    if (nrow(df_gff) == 0) {
        return(list(errmsg=paste0("Error: ", busco, " GFF entry is not found for ", prefix, ". Skipped.")))
    }

    # set up variable
    df_gff_exon <- NULL

    # extract the respective gene index
    gene_idx <- which(df_gff$feature == "gene" & df_gff$start == coordinates$start & df_gff$end == coordinates$stop)

    # extract all exons based on the coordinates
    if (length(gene_idx) != 1) {
        df_gff_exon <- df_gff[df_gff$feature=="exon" & df_gff$start>=coordinates$start & df_gff$end<=coordinates$stop]
    } else {
        # extract all gene indices
        ls_gene_idx <- which(df_gff$feature == "gene")
        entry_idx <- match(gene_idx, ls_gene_idx)

        # subset the GFF table
        df_gff_subset <- NULL
        if (entry_idx == length(ls_gene_idx)) {
            df_gff_subset <- df_gff[gene_idx:nrow(df_gff)]
        } else {
            df_gff_subset <- df_gff[gene_idx:ls_gene_idx[entry_idx+1]-1]
        }

        df_gff_exon <- df_gff_subset[df_gff_subset$feature=="exon"]
    }

    # return error message in no exon found
    if (nrow(df_gff_exon) == 0) {
        return(list(errmsg=paste0("Error: no exon found for ", busco, " for ", prefix, ". Skipped.")))
    }

    # return error message if no intron found
    if (nrow(df_gff_exon) == 1) {
        return(list(errmsg=paste0("Error: no intron found for ", busco, " for ", prefix, ". Skipped.")))
    }

    # remove multiple genes
    df_gff_exon <- f_check_target_id(df_gff_exon)
    # df_gff_exon <- df_gff_exon[order(df_gff_exon$attribute),]

    # add intron entries
    df_output <- data.table::data.table(seqname=character(), source=character(), feature=character(),
                                        start=numeric(), end=numeric(), score=character(), strand=character(),
                                        frame=character(), attribute=character())
    for (i in 1:(nrow(df_gff_exon)-1)) {
        intron_entry <- f_extract_intron_gff(df_gff_exon[i,], df_gff_exon[i+1,], i-1, busco, coordinates$strand)

        # update output data.frame
        df_output <- rbind(df_output, intron_entry)
    }

    # save the new GFF file
    data.table::fwrite(df_output, file=fn_out, sep="\t", quote=F, row.names=F, col.names=F)
}

# function: extract all BUSCO alignments from GFF
f_extract_fasta_from_gff <- function(fn_input, fn_gff, fn_out, fn_prefix, exe_gffread){
    # run gffread
    cmd_gffread <- paste(exe_gffread, "-g", fn_input, "-u", fn_out, fn_gff)
    system(cmd_gffread)

    # open the FASTA sequences
    all_taxa <- Biostrings::readDNAStringSet(fn_out)
    all_headers <- stringr::str_sort(names(all_taxa), numeric=T)

    # save introns as individual FASTA file
    for (i in all_headers) {
        # extract header
        ls_header <- unlist(strsplit(i, split="_"))
        header <- paste0(ls_header[length(ls_header)-1], "_", ls_header[length(ls_header)])

        # output file
        fn_temp_out <- paste0(fn_prefix, "_", header, ".fna")
        
        # save the FASTA in a file
        fasta_seq <- Biostrings::DNAStringSet(all_taxa[[i]])
        names(fasta_seq) <- header
        Biostrings::writeXStringSet(fasta_seq, filepath=fn_temp_out)
    }
}

# function: combine individual FASTA as MSA
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

###################################################################
######                         SETUP                         ######
###################################################################

library(doSNOW)

# create log file
fn_log <- paste0(dir_out, "/log.txt")
log_appender <- log4r::file_appender(fn_log, append = TRUE, layout = log4r::default_log_layout())
fn_logger <- log4r::logger(threshold = "INFO", appenders = log_appender)

unlink(fn_log)

# read input file
df_input <- data.table::fread(file_metadata)

# check if there are multiple refs
ref <- unique(df_input$name[df_input$type=="ref"])
if (length(ref) != 1) {
    log4r::error(fn_logger, "Error: more than one reference sequence in the metadata. Exited.")
    stop()
}

# create output directory
dir_out_ref <- paste0(dir_out, "/", ref, "/")
if (!dir.exists(dir_out_ref)) {
    dir.create(dir_out_ref, recursive=T)
}

dir_out_tree <- paste0(dir_out, "/tree/")
if (!dir.exists(dir_out_tree)) {
    dir.create(dir_out_tree, recursive=T)
}

# output file
fn_robust_tree <- paste0(dir_out, "/robust_tree.tre")
fn_robust_tree_nwk <- paste0(dir_out, "/robust_tree.nwk")
fn_astral_nwk <- paste0(dir_out, "/astral.nwk")

###################################################################
######               EXTRACT REFERENCE BUSCO                 ######
###################################################################

# extract list of reference BUSCOs
busco_ids <- list.files(dir_busco_ref, pattern = "*.fna$", full.names = F, recursive = F)
busco_ids <- sapply(busco_ids, function(x) { gsub(".fna", "", x) })

# extract reference file location
file_fas <- df_input$fasta_dir[df_input$name==ref]
if (!file.exists(file_fas)) {
    log4r::error(fn_logger, paste0("File not found: FASTA file for ", ref, ". Exited."))
    stop()
}

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over shared BUSCOs
foreach (busco = busco_ids) %dopar% {
    require(data.table)
    require(Biostrings)
    require(stringr)

    # BUSCO files
    file_busco_fas <- paste0(dir_busco_ref, "/", busco, ".fna")
    file_busco_gff <- paste0(dir_busco_ref, "/", busco, ".gff")
    if (!file.exists(file_busco_gff)) {
        log4r::warn(fn_logger, paste0("Warn: GFF file for ", busco, " is unavailable. Skipped."))
        return(NULL)
    }
    
    # extract header
    busco_header <- system(paste("grep '^>'", file_busco_fas), intern=T)

    # output files
    file_busco_prefix <- paste0(dir_out_ref, busco)
    file_busco_out <- paste0(file_busco_prefix, ".intron.fas")
    file_busco_gff_out <- paste0(file_busco_prefix, ".intron.gff")

    # extract BUSCO coordinates
    coordinates <- f_extract_coordinates(busco_header, busco, ref)
    if (!is.null(coordinates$errmsg)) {
        f_write_log(fn_log=fn_log, msg=coordinates$errmsg)
        return(NULL)
    }

    # manipulate GFF file
    msg <- f_manipulate_gff(file_busco_gff, coordinates, busco, ref, file_busco_gff_out)
    if (!is.null(msg$errmsg)) {
        f_write_log(fn_log=fn_log, msg=msg$errmsg)
        return(NULL)
    }
    
    # extract BUSCO
    f_extract_fasta_from_gff(file_fas, file_busco_gff_out, file_busco_out, file_busco_prefix, exe_gffread)
}

log4r::info(fn_logger, paste0("File created/modified: BUSCO FASTA alignments for ", ref, "."))
stopCluster(nwcl)

###################################################################
######             EXTRACT MAPPED READS BUSCO                ######
###################################################################

# check mapped reads
ls_prefix <- unique(df_input$name[df_input$type!="ref"])
if (length(ls_prefix) < 1) {
    log4r::error(fn_logger, "Error: less than one mapped reads in the metadata. Exited.")
    stop()
}

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over short reads
for (prefix in ls_prefix) {
    # create outdir
    dir_busco_out <- paste0(dir_out, "/", prefix, "/")
    if (!dir.exists(dir_busco_out)) {
        dir.create(dir_busco_out, recursive=T)
    }

    # check input files
    file_fas <- df_input$fasta_dir[df_input$name==prefix]
    if (!file.exists(file_fas)) {
        log4r::warn(fn_logger, paste0("File not found: FASTA file for ", prefix, ". Skipped."))
        next
    }

    # iterate over shared BUSCOs
    foreach (busco = busco_ids) %dopar% {
        require(Biostrings)
        require(stringr)

        # input files
        file_busco_gff_out <- paste0(dir_out_ref, "/", busco, ".intron.gff")
        if (!file.exists(file_busco_gff_out)) {
            log4r::warn(fn_logger, msg=paste0("File not found: ", busco, " GFF for ", ref, ". Skipped."))
            return(NULL)
        }

        # output files
        file_busco_prefix <- paste0(dir_busco_out, busco)
        file_busco_out <- paste0(file_busco_prefix, ".intron.fas")

        # extract BUSCO
        f_extract_fasta_from_gff(file_fas, file_busco_gff_out, file_busco_out, file_busco_prefix, exe_gffread)
    }

    log4r::info(fn_logger, paste0("File created/modified: BUSCO FASTA alignments for ", prefix, "."))
}

stopCluster(nwcl)

###################################################################
######                         MSA                           ######
###################################################################

# extract list of reference BUSCO introns
intron_ids <- list.files(dir_out_ref, pattern = "*.fna$", full.names = F, recursive = F)
intron_ids <- sapply(intron_ids, function(x) { gsub(".fna", "", x) })

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over shared BUSCOs
foreach (busco = intron_ids) %dopar% {
    dir_output <- paste0(dir_out_tree, busco, "/")
    if (!dir.exists(dir_output)) {
        dir.create(dir_output, recursive=T)
    }

    # output MSA file
    fn_out <- paste0(dir_output, busco, ".fas")

    # extract BUSCO alignment from reference sequences
    file_busco_ref <- paste0(dir_out_ref, busco, ".fna")
    if (!file.exists(file_busco_ref)) {
        next
    }

    f_fasta2msa(file_busco_ref, ref, fn_out)

    # extract BUSCO alignment from short reads
    for (prefix in ls_prefix) {
        file_busco_read <- paste0(dir_out, "/", prefix, "/", busco, ".fna")
        if (!file.exists(file_busco_read)) {
            next
        }

        f_fasta2msa(file_busco_read, prefix, fn_out)
    }
}

stopCluster(nwcl)

###################################################################
######                     INTRON TREE                       ######
###################################################################

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over shared BUSCOs
ls_trees <- foreach (busco = intron_ids, .combine='c') %dopar% {
    dir_output <- paste0(dir_out_tree, busco, "/")

    # output IQ-Tree2 files
    fn_fasta <- paste0(dir_output, busco, ".fas")
    if (file.exists(fn_fasta)) {
        # run IQ-Tree2
        cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_fasta,
                         "-bb 1000",
                         "-T 1 --quiet -redo")
        system(cmd_iqtree2)
    }

    # extract only tree with >95 bootstrap value
    fn_treefile <- paste0(dir_output, busco, ".fas.treefile")
    if (file.exists(fn_treefile)) {
        tl <- ape::read.tree(fn_treefile)

        # extract bootstrap values
        bl <- subset(tl$node.label, tl$node.label != "")
        if (!is.null(bl) && mean(as.numeric(bl)) >= 95) {
            return(list(busco=busco, tree=ape::write.tree(tl)))
        }
    }
}

# save the output in a file
f_write_log(fn_robust_tree, msg=unlist(ls_trees[names(ls_trees)=="busco"]))
f_write_log(fn_robust_tree_nwk, msg=unlist(ls_trees[names(ls_trees)=="tree"]))

log4r::info(fn_logger, "File created/modified: BUSCO intron trees.")
stopCluster(nwcl)

###################################################################
######                      ASTRAL-III                       ######
###################################################################
cmd_astral <- paste("java -jar", exe_astral,
                    "-i", fn_robust_tree_nwk,
                    "-o", fn_astral_nwk,
                    "-t 2")
system(cmd_astral)