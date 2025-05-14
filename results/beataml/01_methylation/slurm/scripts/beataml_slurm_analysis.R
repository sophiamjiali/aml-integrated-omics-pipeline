library(doParallel)
library(minfi)
library(doRNG)

# fetch home directories from the environment variables
home_dir <- Sys.getenv("HOME")
scratch_dir <- Sys.getenv("SCRATCH")

# initialize clusters for parallel computation
cl <- makeCluster(40)
registerDoParallel(cl)

# fetch the generated input data
data <- readRDS(file.path(home_dir, "slurm", "data", "beataml_slurm_data.rds"))
info <- readRDS(file.path(home_dir, "slurm", "data", "beataml_slurm_info.rds"))
design <- readRDS(file.path(home_dir, "slurm", "data", "beataml_slurm_design.rds"))

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
  pickCutoffQ = 0.99,
  maxGap = 500,
  minNum = 3,
  nullMethod = "permutation",
  B = 1000,
  smooth = TRUE,
  smoothFunction = loessByCluster,
  useWeights = FALSE,
  verbose = TRUE)

# save the DMRs to the SCRATCH output directory
saveRDS(dmrs, file.path(scratch_dir, "beataml_raw_dmrs.rds"))
stopCluster(cl)
