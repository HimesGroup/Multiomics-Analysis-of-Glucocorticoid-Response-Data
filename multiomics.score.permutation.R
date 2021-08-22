# /usr/bin/R
# author: Mengyuan Kan (mengykan@pennmedicine.upenn.edu)
# Licensed under the MIT License

library(argparse)

###
# Arguments
###

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--score_fn", type="character",
                    help="Filename of multiomics integrative scores with full path.")
parser$add_argument("--out_fn", type="character",
                    help="Output filename of scores for GR-binding.")
parser$add_argument("--n_core", type="integer", default=20,
                    help="Number of threads used to run permutation [default %(default)s]")
parser$add_argument("--nperm", type="integer", default=1000000,
                    help="Number of permutations [default %(default)s]")

args <- parser$parse_args()



###
# Check Arguements
###

score_fn <- args$score_fn
out_fn <- args$out_fn
n_core <- args$n_core
nperm <- args$nperm


#score_fn <- "/project/bhimeslab//GSK_share/multiomics/GSK.score.txt"
#n_core = 20
#nperm = 1000000
#out_fn <- "/project/bhimeslab/GSK_share/multiomics/GSK.score.perm.res"

# required arguements
if (is.null(score_fn)|is.null(out_fn)) {
  cat("Some required arguements are missing. Please check!\n")
  parser$print_help()
  quit()
}
# check if file exists
if (!file.exists(score_fn)) {
  cat("Some pre-defined files do not exist.\n")
  quit()
}

out_path <- dirname(out_fn)
if (!dir.exists(out_path)) {
  cat("Output directory", out_path, "does not exist.\n")
  quit()
}

# load libraries
library(dplyr)
library(parallel)

dat.score <- read.table(score_fn, header=T, sep="\t")
rownames(dat.score) <- dat.score$SNP
dat.score$SNP <- NULL
score.org <- apply(dat.score, 1, sum)
# permute dataframe function
permvec_func <- function(vec){vec[!is.na(vec)] <- sample(vec[!is.na(vec)]); return(vec)} 
perm_func <- function(rand){
  set.seed(rand)
  dat.tmp=apply(dat.score,2,permvec_func)
  score_perm=apply(dat.tmp, 1, sum)
  compare = as.numeric(score_perm >= score.org)
  return(compare)
}
rands = c(1:nperm)
list.perm <- mclapply(rands, perm_func, mc.cores = n_core)
perm.count = Reduce(`+`, list.perm)
perm.pval <- (perm.count+1)/(nperm+1)
perm.qval <- p.adjust(perm.pval, method = "BH")
res <- data.frame(SNP=names(score.org), score.org, perm.count, perm.pval, perm.qval) %>% dplyr::arrange(perm.pval)
write.table(res, out_fn, row.names=F, col.names=T, sep="\t", quote=F)

print(sessionInfo())

