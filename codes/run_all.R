#################################
codedir <- "~/BusIER/codes"
prefix <- ""
outdir <- "~/busco"
thread <- 5
redo <- TRUE

# data download
file_refseq <- ""
file_shortreads <- ""

exe_datasets <- "datasets"
dir_sratoolkit <- "sratoolkit/bin/"

# BUSCO check
exe_bwamem2 <- "bwa-mem2"
exe_samtools <- "samtools"
exe_bcftools <- "bcftools"

exe_busco <- "busco"
exe_iqtree2 <- "iqtree2"
exe_mafft <- "mafft"

n_short_reads <- 3

busco_lineage <- "eudicots_odb10"
busco_mode <- "genome"
#################################

# set up outdir
outdir_prefix <- paste0(outdir, "/", prefix, "/")
if (!dir.exists(outdir_prefix)) {
    dir.create(outdir_prefix, recursive=T)
}

# data download
rmarkdown::render(input=paste0(codedir,"/1_data_download/1_main.Rmd"),
                  output_file=paste0(outdir_prefix, prefix, ".download.html"),
                  params=list(prefix=prefix, codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                              file_refseq=file_refseq, file_shortreads=file_shortreads,
                              exe_datasets=exe_datasets, dir_sratoolkit=dir_sratoolkit),
                  quiet=TRUE)

# BUSCO check
rmarkdown::render(input=paste0(codedir,"/2_busco_check/1_main.Rmd"),
                  output_file=paste0(outdir_prefix, prefix, ".check.html"),
                  params=list(prefix=prefix, codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                              exe_bwamem2=exe_bwamem2, exe_samtools=exe_samtools, exe_bcftools=exe_bcftools, exe_busco=exe_busco, exe_iqtree2=exe_iqtree2, exe_mafft=exe_mafft,
                              n_short_reads=n_short_reads, busco_lineage=busco_lineage, busco_mode=busco_mode),
                  quiet=TRUE)

#################################