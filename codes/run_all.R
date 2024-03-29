#################################
codedir <- "~/BusIER/codes"
prefix <- ""
outdir <- "~/busco"
thread <- 5
redo <- TRUE

# data download
file_refseq <- paste0(codedir, "/../data/eucs_refseq.txt")
file_shortreads <- paste0(codedir, "/../data/eucs_shortreads.txt")
file_adapters <- paste0(codedir, "/../data/eucs_adapter.txt")
file_refseq_chr <- paste0(codedir, "/../data/eucs_refseq_chr.txt")

exe_datasets <- "datasets"
bin_sratoolkit <- "sratoolkit/bin/"

# BUSCO check
exe_adapterremoval <- "AdapterRemoval"
exe_bwamem2 <- "bwa-mem2"
exe_samtools <- "samtools"
exe_bcftools <- "bcftools"

exe_gffread <- "gffread"

exe_busco <- "busco"
exe_iqtree2 <- "iqtree2"
exe_mafft <- "mafft"

busco_lineage <- "eudicots_odb10"
busco_mode <- "genome"
type <- "coordinate"
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
                              exe_datasets=exe_datasets, bin_sratoolkit=bin_sratoolkit),
                  quiet=TRUE)

# BUSCO check
rmarkdown::render(input=paste0(codedir,"/2_busco_check/1_main.Rmd"),
                  output_file=paste0(outdir_prefix, prefix, ".check.html"),
                  params=list(prefix=prefix, codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                              exe_adapterremoval=exe_adapterremoval, exe_bwamem2=exe_bwamem2, exe_samtools=exe_samtools, exe_bcftools=exe_bcftools,
                              exe_gffread=exe_gffread, exe_busco=exe_busco, exe_iqtree2=exe_iqtree2, exe_mafft=exe_mafft,
                              file_refseq=file_refseq, file_shortreads=file_shortreads, file_adapters=file_adapters, file_refseq_chr=file_refseq_chr,
                              busco_lineage=busco_lineage, busco_mode=busco_mode, type=type),
                  quiet=TRUE)

#################################