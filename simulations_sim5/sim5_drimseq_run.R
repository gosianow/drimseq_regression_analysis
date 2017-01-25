##############################################################################

# BioC 3.4
# Created 23 Jan 2017 

##############################################################################
Sys.time()
##############################################################################

library(BiocParallel)
library(DRIMSeq)
library(limma)

# devtools::load_all("/Users/gosia/Dropbox/UZH/package_devel/DRIMSeq")

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/drimseq_regression_analysis/simulations_sim5/drosophila_node_nonull'
out_dir='drimseq_1_3_3/kallistoprefiltered5'
workers=3
count_method='kallistoprefiltered5'

one_way=FALSE
coef_mode='optim'
disp_moderation='trended'
out_prefix='drimseq_reg_optim_trended_'

one_way=TRUE
coef_mode='optim'
disp_moderation='trended'
out_prefix='drimseq_oneway_trended_'


##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

cat(paste0(args, collapse = "\n"), fill = TRUE)

##############################################################################

setwd(rwd)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}

########################################################
# Create metadata file
########################################################

metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata


##########################################################################
# Load counts
##########################################################################


if(count_method == "htseq")
  count_dir <- "2_counts/dexseq_nomerge/dexseq"
if(count_method == "kallisto")
  count_dir <- "2_counts/kallisto/kallisto"
if(count_method == "htseqprefiltered15")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast15/dexseq"
if(count_method == "htseqprefiltered5")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast5/dexseq"
if(count_method == "kallistofiltered5")
  count_dir <- "2_counts/kallisto_txfilt_5/kallisto"
if(count_method == "kallistoprefiltered5")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/kallisto_kallistoest_atleast5/kallisto"

count_dir

### load counts
counts_list <- lapply(1:6, function(i){
  # i = 1
  cts <- read.table(paste0(count_dir, i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(cts) <- c("group_id", paste0("sample_", i))  
  return(cts)
})

counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)

counts <- counts[!grepl(pattern = "_", counts$group_id),]

counts[, c("gene_id", "feature_id")] <- strsplit2(counts[,1], ":")

counts <- counts[, c("gene_id", "feature_id", metadata$sample_id)]

head(counts)

##########################################################################
# DRIMSeq analysis
##########################################################################


# Create the dmDSdata object
d <- dmDSdata(counts = counts, samples = metadata)
d


# Filtering
d <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0)

# Create the design matrix
design <- model.matrix( ~ group, data = samples(d))


# Calculate dispersion
set.seed(1234)
d <- dmDispersion(d, design = design, one_way = one_way, coef_mode = coef_mode, disp_moderation = disp_moderation, BPPARAM = BPPARAM)


ggp <- plotDispersion(d)

pdf(file.path(out_dir, paste0(out_prefix, "dispersion.pdf")))
print(ggp)
dev.off()


# Fit full model proportions
d <- dmFit(d, design = design, one_way = one_way, bb_model = one_way, coef_mode = coef_mode, BPPARAM = BPPARAM)


# Fit null model proportions and test for DS
d <- dmTest(d, coef = 2, one_way = one_way, bb_model = one_way, coef_mode = coef_mode, BPPARAM = BPPARAM)


res_gene <- results(d, level = "gene")

write.table(res_gene, file.path(out_dir, paste0(out_prefix, "results_gene.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



if(one_way){
  
  res_feature <- results(d, level = "feature")
  
  res_twostage01 <- dmTwoStageTest(d, FDR = 0.01, BPPARAM = BPPARAM)
  res_twostage05 <- dmTwoStageTest(d, FDR = 0.05, BPPARAM = BPPARAM)
  res_twostage10 <- dmTwoStageTest(d, FDR = 0.10, BPPARAM = BPPARAM)
  
  write.table(res_feature, file.path(out_dir, paste0(out_prefix, "results_feature.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  write.table(res_twostage01, file.path(out_dir, paste0(out_prefix, "results_twostage01.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(res_twostage05, file.path(out_dir, paste0(out_prefix, "results_twostage05.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(res_twostage10, file.path(out_dir, paste0(out_prefix, "results_twostage10.txt")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
}



save(d, file = file.path(out_dir, paste0(out_prefix, "d.RData")))




sessionInfo()







