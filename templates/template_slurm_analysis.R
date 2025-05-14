library(doParallel)
library(minfi)
library(doRNG)

# fetch home directories from the environment variables
home_dir <- Sys.getenv("HOME")
scratch_dir <- Sys.getenv("SCRATCH")

# initialize clusters for parallel computation
cl <- makeCluster({{n_cores}})
registerDoParallel(cl)

# fetch the generated input data
data <- readRDS(file.path(home_dir, "slurm", "data", "{{data_file}}"))
info <- readRDS(file.path(home_dir, "slurm", "data", "{{info_file}}"))
design <- readRDS(file.path(home_dir, "slurm", "data", "{{design_file}}"))

# generate DMRs via bumphunter
set.seed(123)
RNGkind("L'Ecuyer-CMRG")
options(rngtools.RNGseed = 123)

dmrs <- bumphunter(
  object = data,
  design = design,
  chr = info$chr,
  pos = info$pos,
  coef = 2,
  cutoff = NULL,
  pickCutoff = TRUE,
  pickCutoffQ = {{cutoff_q}},
  maxGap = {{max_gap}},
  minNum = {{min_probes}},
  nullMethod = "permutation",
  B = {{permutations}},
  smooth = TRUE,
  smoothFunction = loessByCluster,
  useWeights = FALSE,
  verbose = TRUE)

# save the DMRs to the SCRATCH output directory
saveRDS(dmrs, file.path(scratch_dir, "{{dmrs_file}}"))
stopCluster(cl)
