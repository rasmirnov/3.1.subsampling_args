configfile: "config.yaml"

# link на исходные скрипты с scclustevel: https://github.com/crazyhottommy/pyflow_seurat_parameter

NUM_OF_SUBSAMPLE = config["num_of_subsample"]
rate = config["subsample_rate"]
resolutions = config["subsample_resolutions"].strip().split()
INPUT_SEURAT = config["input_seurat"]
rscript_path = config['rscript']


SUBSAMPLE_RESOLUTION = expand("subsample/subsample_resolution_{resolution}_round_{run_id}.rds", resolution = resolutions, run_id = range(NUM_OF_SUBSAMPLE))

FULLSAMPLE_RESOLUTION = expand("full_sample_preprocess/full_sample_resolution_{resolution}.rds", resolution = resolutions)

TARGETS = []

TARGETS.extend(SUBSAMPLE_RESOLUTION)
TARGETS.append("gather_subsample.rds")
TARGETS.append("gather_full_sample.rds")

localrules: all, gather_subsample, gather_full_sample_preprocess
rule all:
    input: TARGETS

## the full data set, preprocessing using a set of k, resolution and PC
       
rule full_sample_preprocess:
        input: rda=INPUT_SEURAT
        output: rds="full_sample_preprocess/full_sample_resolution_{resolution}.rds"
        log: "00log/full_sample_resolution_{resolution}.log"
        shell: "{rscript_path} scripts/preprocess.R --data {input.rda} --res {wildcards.resolution}"


rule gather_full_sample_preprocess:
        input: rds = FULLSAMPLE_RESOLUTION
        output: rds="gather_full_sample.rds"
        log: "00log/full_sample_gather_idents.log"
        shell: "{rscript_path}  scripts/gather_fullsample.R --data {input.rds} --out_file {output.rds}"


## subsample e.g. 80% of the cells and re-do the clustering for n times
rule subsample_cluster:
        input: rds="full_sample_preprocess/full_sample_resolution_{resolution}.rds"
        output: "subsample/subsample_resolution_{resolution}_round_{run_id}.rds"
        log: "00log/subsample_resolution_{resolution}_round_{run_id}.log"
        shell: """
               {rscript_path}  scripts/subsample.R \
               --data {input.rds} \
               --res {wildcards.resolution} \
               --run_id {wildcards.run_id} \
               --rate {rate} 2> {log}\
               """


## gather the subsampled and reclustered cell idents
rule gather_subsample:
        input: rds = SUBSAMPLE_RESOLUTION
        output: "gather_subsample.rds"
        log: "00log/gather_subsample.log"
        shell: "{rscript_path}  scripts/gather_subsample.R --data {input.rds}"
