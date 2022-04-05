library(tidyverse)
library(argparse)
library(Seurat)

## SET VARIABLES

parser <-  ArgumentParser()
parser$add_argument('--data',
                    type = "character",
                    help = 'path to rdata object')

args <- parser$parse_args()

print(args)
rdss<- args$data

get_df<- function(rds){
	res<- readRDS(rds)
	return(res)
}

dat.list<- lapply(rdss, get_df)

gather_idents<- do.call(bind_rows, dat.list)
saveRDS(gather_idents, file = "gather_subsample.rds")
