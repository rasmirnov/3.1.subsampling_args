library(Seurat)
library(tidyverse)
library(argparse)

## SET VARIABLES

parser <-  ArgumentParser()
parser$add_argument('--data',
                    type = "character",
                    help = 'path to rdata object')
parser$add_argument('--out_file',
                    type = "character",
                    help = 'output filename')

args <- parser$parse_args()

print(args)

seurat_obj<- get(load(args$data))
rdss<- args$data

get_idents<- function(rds){
        x<- readRDS(rds)
        resolution<- gsub("full_sample_resolution_([0-9\\.]+).rds", "\\1", basename(rds))
        df<- tibble::tibble(resolution = resolution, original_ident_full = list(Idents(x)))
        return(df)
}

dat.list<- lapply(rdss, get_idents)

gather_idents<- do.call(bind_rows, dat.list)
saveRDS(gather_idents, file = args$out_file)
