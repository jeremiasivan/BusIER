#################################
codedir <- "~/BusIER/codes"
prefix <- ""
outdir <- "~/busco"
thread <- 5
redo <- TRUE

# data download
file_refseq <- ""
file_shortreads <- paste0(codedir, "/../data/eucs_shortreads.txt")
file_adapters <- paste0(codedir, "/../data/eucs_adapter.txt")

exe_datasets <- "datasets"
bin_sratoolkit <- "sratoolkit/bin/"

# gene tree
gene_type <- "easy353"

exe_adapterremoval <- "AdapterRemoval"
exe_build_database <- "build_database.py"
exe_easy353 <- "easy353.py"

easy353_taxonomy <- "Myrtaceae"
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
                              file_refseq=file_refseq, file_shortreads=file_shortreads, file_adapters=file_adapters,
                              exe_datasets=exe_datasets, bin_sratoolkit=bin_sratoolkit),
                  quiet=TRUE)

# BUSCO check
rmarkdown::render(input=paste0(codedir,"/3_gene_tree_run/1_main.Rmd"),
                  output_file=paste0(outdir_prefix, prefix, ".genetree.html"),
                  params=list(prefix=prefix, codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                              exe_adapterremoval=exe_adapterremoval, exe_build_database=exe_build_database, exe_easy353=exe_easy353,
                              gene_type=gene_type, easy353_taxonomy=easy353_taxonomy),
                  quiet=TRUE)

#################################