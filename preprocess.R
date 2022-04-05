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

args <- parser$parse_args()

print(args)

seurat_obj<- get(load(args$data))
resolution<- args$res


PreprocessSubsetData<- function(object,
                                variable.features.n = 3000,
                                num.pc = 20,
                                pc.use = NULL,
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
        
        ## use future for parallelization
        ##future::plan("multiprocess", workers = workers)
        #
        meta.data.colnames<- object@meta.data %>% colnames()
        vars.to.regress<- c("percent.mt","nFeature_RNA")
        # in case the seurat object does not have percent.mito in metadata
        vars.to.regress<- vars.to.regress[vars.to.regress %in% meta.data.colnames]
        # default is on variable features only, omit the features argument
        # SCTransform replaces NormalizeData, ScaleData and FindVariableFeatures
        
        if(useSCTransform==TRUE){
                object<- SCTransform(object, vars.to.regress = vars.to.regress,
                                     variable.features.n = variable.features.n, verbose = FALSE)  
        }else{
                stop("The SCTransform method for normalization is the only method supported by this function.  If you wish to use the approach that involves NormalizeData, ScaleData, and FindVariableFeatures and enables use of the Jackstraw procedure for determining which PCs to use please use the PreprocessSubsetDataV2 function from the scclusteval R package.")
                
        }

        object<- RunPCA(object = object, features = VariableFeatures(object = object),
                        npcs = num.pc)

        
        if(is.null(pc.use)){
                pc.use <- num.pc
                message("SCTransform is being used and the Jackstraw procedure for determining which PCs to use is not compatable with this procedure. Since pc.use was not specified it is being automatically set to num.pc")
        }

        # add significant pc number to metadata, need to have names same as the cells
        pc.use.meta<- rep(pc.use, length(colnames(object)))
        names(pc.use.meta)<- colnames(object)
        object<- AddMetaData(object = object, metadata = pc.use.meta, col.name = "pc.use")
        object<- FindNeighbors(object, dims = 1:pc.use, k.param = k.param, nn.eps = nn.eps,
                               verbose = FALSE, reduction = "pca", force.recalc = TRUE)
        object <- FindClusters(object = object,
                                n.start = n.start,
                                resolution = resolution,
                                verbose = FALSE)
        return(object)
}


# PreprocessSubsetData_pars<- snakemake@params[["PreprocessSubsetData_pars"]]
## this is not subsetted data, but the PreprocessSubsetData function can be used as well for any seurat object
# 'eval' evaluetes the specific expression and if it is a Python readable --> execute it
#? why parse? --> converts smt to a dataframe?
seurat_obj<- eval(parse(text=paste("PreprocessSubsetData", "(", "seurat_obj,",
                                   "resolution=", resolution, ")")))
saveRDS(seurat_obj, file = paste0("full_sample_preprocess/full_sample_", "resolution_", resolution, ".rds"))
