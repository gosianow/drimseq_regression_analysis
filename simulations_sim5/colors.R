##############################################################################

# BioC 3.4
# Created 23 Jan 2017 

##############################################################################

Sys.time()

##############################################################################

library(ggplot2)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/drimseq_regression_analysis/simulations_sim5'
out_dir='/Users/gosia/Dropbox/UZH/drimseq_regression_analysis/simulations_sim5'


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

dir.create(out_dir, recursive = TRUE)

##############################################################################

methods <- c("dexseq", "drimseq_reg_optim_none", "drimseq_reg_optim_trended", "drimseq_oneway_none","drimseq_oneway_trended")

methods <- factor(methods, levels = methods)


colors_df <- data.frame(methods = methods, colors = c("royalblue3", "mediumvioletred", "orchid4", "orange2", "chocolate4"), colors_hex = c("#3A5FCD", "#C71585", "#8B4789", "#EE9A00", "#8B4513"))
colors_df$colors <- as.character(colors_df$colors)
colors_df$colors_hex <- as.character(colors_df$colors_hex)


ggp <- ggplot(colors_df, aes(x = methods, y=rep(1, length(colors)), fill = methods)) + geom_bar(stat="identity") + scale_fill_manual(values = colors_df$colors) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


pdf(paste0(out_dir, "/colors.pdf"), 10, 5)
print(ggp)
dev.off()


colors <- colors_df$colors
names(colors) <- colors_df$methods



write.table(colors_df, paste0(out_dir, "/colors.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(colors, colors_df, file = paste0(out_dir, "/colors.Rdata"))








