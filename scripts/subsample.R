library(Seurat)
library(tidyverse)
library(argparse)

## SET VARIABLES

parser <-  ArgumentParser()
parser$add_argument('--data',
                    type = "character",
                    help = 'path to rdata object')
parser$add_argument('--res',
                    type = "character",
                    help = 'resolution')
parser$add_argument('--run_id',
                    type = "character",
                    help = 'run_id')
parser$add_argument('--rate',
                    type = "character",
                    help = 'rate')

args <- parser$parse_args()

print(args)

seurat_obj<- get(load(args$data))
resolution<- args$res
rdss<- args$data
run_id<- args$run_id

RandomSubsetData<- function(object, rate, random.subset.seed = NULL, ...){
        ncells<- nrow(object@meta.data)
        ncells.subsample<- round(ncells * rate)

        set.seed(random.subset.seed)

        selected.cells<- sample(colnames(object), ncells.subsample)
        object<- subset(object, cells =  selected.cells,
                            ...)
        return(object)
}

subset_seurat_obj<- RandomSubsetData(seurat_obj, rate = args$rate)
original_ident<- Idents(subset_seurat_obj)

PreprocessSubsetData<- function(object,
                                variable.features.n = 3000,
                                num.pc = 20,
                                pc.use = 50,
                                n.start = 100,
                                nn.eps = 0,
                                resolution = 0.8,
                                k.param = 30,
                                useSCTransform = TRUE,
                                ...){
    if(!is.null(pc.use)){
        if(pc.use > num.pc){
            stop("Specify the maximum pc.use number as less than or equal to the total num.pc")
        }
    }

        meta.data.colnames<- object@meta.data %>% colnames()
        vars.to.regress<- c("percent.mt","nFeature_RNA")
        vars.to.regress<- vars.to.regress[vars.to.regress %in% meta.data.colnames]
        if(useSCTransform==TRUE){
            object<- SCTransform(object, vars.to.regress = vars.to.regress,
                                 variable.features.n = variable.features.n, verbose = FALSE)  
        }else{
            stop("The SCTransform method for normalization is the only method currently supported by this function.  If you wish to use the approach that involves NormalizeData, ScaleData, and FindVariableFeatures and enables use of the Jackstraw procedure for determining which PCs to use please use the PreprocessSubsetDataV2 function from the scclusteval R package.")
            
        }

        object<- RunPCA(object = object, features = VariableFeatures(object = object),
                        npcs = num.pc)

  
        if(is.null(pc.use)){
            pc.use <- num.pc
            message("SCTransform is being used and the Jackstraw procedure for determining which PCs to use is not compatable with this procedure. Since pc.use was not specified it is being automatically set to num.pc")
        }
     pc.use.meta<- rep(pc.use, length(colnames(object)))
        names(pc.use.meta)<- colnames(object)
        object<- AddMetaData(object = object, metadata = pc.use.meta, col.name = "pc.use")
        object<- FindNeighbors(object, dims = 1:pc.use,
                               verbose = FALSE, reduction = "pca", force.recalc = TRUE)
        object <- FindClusters(object = object,
                                resolution = resolution,
                                verbose = FALSE)
        return(object)
}

## after reprocessing, the ident slot will be updated with the new cluster id
command<- paste("PreprocessSubsetData", "(", "subset_seurat_obj,", 
                                   "resolution=", resolution, ")")
subset_seurat_obj<- eval(parse(text=command))

res<- tibble::tibble(resolution = resolution,original_ident = list(original_ident),
    recluster_ident = list(Idents(subset_seurat_obj)), round = run_id)


outfile<- paste0("subsample/subsample_", "_resolution_", resolution, "_round_", run_id, ".rds")
saveRDS(res, file = outfile)

## make sure it is not empty file
info<- file.info(outfile)
if (info$size == 0) {
    quit(status = 1)
}
