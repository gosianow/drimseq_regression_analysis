##############################################################################

# BioC 3.4
# Created 23 Jan 2017

##############################################################################
Sys.time()
##############################################################################

library(iCOBRA)
library(Hmisc)
library(DEXSeq)
library(DRIMSeq)
library(plyr)
library(limma)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/drimseq_regression_analysis/simulations_sim5/drosophila_node_nonull'
out_dir='drimseq_1_3_3_comparison/kallistoprefiltered5'
workers=4
count_method='kallistoprefiltered5'

colors_path='/Users/gosia/Dropbox/UZH/drimseq_regression_analysis/simulations_sim5/colors.Rdata'
drimseq_res_dir='drimseq_1_3_3/kallistoprefiltered5'

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

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# Merge results for iCOBRA
# -----------------------------------------------------------------------------

# Results from DEXSeq

if(count_method == "htseq")
  results_dir <- "4_results/dexseq_htseq_nomerge"
if(count_method == "kallisto")
  results_dir <- "4_results/dexseq_kallisto"
if(count_method == "htseqprefiltered15")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast15"
if(count_method == "htseqprefiltered5")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5"
if(count_method == "kallistofiltered5")
  results_dir <- "4_results/dexseq_kallisto_txfilt_5"
if(count_method == "kallistoprefiltered5")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_kallisto_kallistoest_atleast5"

results_dir


results_padj <- list()


rt <- read.table(paste0(results_dir, ".txt"), header = TRUE, as.is = TRUE)
head(rt)

colnames(rt) <- c("gene_id", "dexseq")
head(rt)

results_padj[["dexseq"]] <- rt


# Results from DRIMSeq

files <- list.files(path = drimseq_res_dir, pattern = "_results_gene.txt" )
files

for(i in 1:length(files)){
  # i = 1
  method_name <- gsub(pattern = "_results_gene.txt", replacement = "", x = files[i])
  
  rt <- read.table(file.path(drimseq_res_dir, files[i]), header = TRUE, as.is = TRUE)
  head(rt)
  
  rt <- rt[,c("gene_id","adj_pvalue")]
  colnames(rt) <- c("gene_id", method_name)
  
  results_padj[[method_name]] <- rt 
  
}


results_padj <- Reduce(function(...) merge(..., by = "gene_id", all=TRUE, sort = FALSE), results_padj)
rownames(results_padj) <- results_padj$gene_id



# Adjust colors

load(colors_path)

colors_df$methods <- as.character(colors_df$methods)

results_padj <- results_padj[, colnames(results_padj) %in% colors_df$methods]

keep_methods <- colors_df$methods %in% colnames(results_padj)

clolors <- colors[keep_methods]
colors_df <- colors_df[keep_methods, , drop = FALSE]

results_padj <- results_padj[, colors_df$methods]



# -----------------------------------------------------------------------------
# Load simulation info
# -----------------------------------------------------------------------------


truth_file <- list.files("3_truth/", pattern = "truth")
truth_file

truth <- read.table(paste0("3_truth/", truth_file), header = TRUE, as.is = TRUE, sep = "\t")

truth <- truth[, c("gene", "ds_status", "de_status", "TPM", "nbr_isoforms", "diff_IsoPct", "nbrexonbins")]
rownames(truth) <- truth$gene


# -----------------------------------------------------------------------------
# Plots with iCOBRA
# -----------------------------------------------------------------------------

cobradata <- COBRAData(padj = results_padj, truth = truth)


### FDR TPR overall

cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", aspects = c("fdrtpr", "fdrtprcurve"), onlyshared = FALSE, maxsplit = Inf)

cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = factor(colors_df$methods, levels = colors_df$methods), colorscheme = colors[basemethods(cobraperf)])

cobraplot <- iCOBRA:::reorder_levels(cobraplot, colors_df$methods)


## Plot the points and the curve

xaxisrange <- c(0, 0.6)
yaxisrange <- c(0.4, 1)

ggp <- plot_fdrtprcurve(cobraplot, plottype = c("curve", "points"), pointsize = 3, xaxisrange = c(0, xaxisrange[2]), yaxisrange = yaxisrange)

ggp <- ggp + 
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 10), strip.text = element_text(size = 11)) + 
  guides(colour = guide_legend(nrow = 2))

pdf(file.path(out_dir, "fdrtprcurve.pdf"), 6, 7)
print(ggp)
dev.off()






# -----------------------------------------------------------------------------
# Feature-level comparison 
# -----------------------------------------------------------------------------


results_feature <- list()

# DEXSeq results

load(paste0(results_dir, ".Rdata"))

results_feature[["dexseq"]] <- data.frame(feature_id = substring(res$featureID, 2, nchar(res$featureID)), dexseq = res$padj, stringsAsFactors = FALSE)


# Plot a histogram of p-values

ggp <- DRIMSeq:::dm_plotPValues(pvalues = res$pvalue) +
  geom_histogram(breaks = seq(0, 1, by = 0.01), fill = "grey50")

pdf(file.path(out_dir, "dexseq_feature_hist_pvalues.pdf"))
print(ggp)
dev.off()



# Results from DRIMSeq

files <- list.files(path = drimseq_res_dir, pattern = "_results_feature.txt" )
files

for(i in 1:length(files)){
  # i = 1
  method_name <- gsub(pattern = "_results_feature.txt", replacement = "", x = files[i])
  
  rt <- read.table(file.path(drimseq_res_dir, files[i]), header = TRUE, as.is = TRUE)
  head(rt)
  
  # Plot a histogram of p-values
  
  ggp <- DRIMSeq:::dm_plotPValues(pvalues = rt$pvalue) +
    geom_histogram(breaks = seq(0, 1, by = 0.01), fill = "grey50")
  
  pdf(file.path(out_dir, paste0(method_name, "_feature_hist_pvalues.pdf")))
  print(ggp)
  dev.off()
  
  
  rt <- rt[,c("feature_id","adj_pvalue")]
  colnames(rt) <- c("feature_id", method_name)
  
  results_feature[[method_name]] <- rt 
  
  
}


results_feature <- Reduce(function(...) merge(..., by = "feature_id", all=TRUE, sort = FALSE), results_feature)
rownames(results_feature) <- results_feature$feature_id



# Adjust colors

load(colors_path)

colors_df$methods <- as.character(colors_df$methods)

results_feature <- results_feature[, colnames(results_feature) %in% colors_df$methods]

keep_methods <- colors_df$methods %in% colnames(results_feature)

clolors <- colors[keep_methods]
colors_df <- colors_df[keep_methods, , drop = FALSE]

results_feature <- results_feature[, colors_df$methods]



# -----------------------------------------------------------------------------
# Load simulation info
# -----------------------------------------------------------------------------

simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE, sep = "\t")

truth_features <- data.frame(ds_status = simulation_details$transcript_ds_status)
rownames(truth_features) <- simulation_details$transcript_id



# -----------------------------------------------------------------------------
# Plots with iCOBRA
# -----------------------------------------------------------------------------

cobradata <- COBRAData(padj = results_feature, truth = truth_features)


### FDR TPR overall

cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", aspects = c("fdrtpr", "fdrtprcurve"), onlyshared = FALSE, maxsplit = Inf)

cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = factor(colors_df$methods, levels = colors_df$methods), colorscheme = colors[basemethods(cobraperf)])

cobraplot <- iCOBRA:::reorder_levels(cobraplot, colors_df$methods)


## Plot the points and the curve

xaxisrange <- c(0, 0.6)
yaxisrange <- c(0.4, 1)

ggp <- plot_fdrtprcurve(cobraplot, plottype = c("curve", "points"), pointsize = 3, xaxisrange = c(0, xaxisrange[2]), yaxisrange = yaxisrange)

ggp <- ggp + 
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 10), strip.text = element_text(size = 11)) + 
  guides(colour = guide_legend(nrow = 2))

pdf(file.path(out_dir, "fdrtprcurve_feature.pdf"), 6, 7)
print(ggp)
dev.off()

























