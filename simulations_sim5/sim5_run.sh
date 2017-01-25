#!/bin/bash
## Define paths to software and reference files

RCODE=/Users/gosia/Dropbox/UZH/drimseq_regression_code/simulations_sim5
RWD=/Users/gosia/Dropbox/UZH/drimseq_regression_analysis/simulations_sim5
ROUT=$RWD/Rout

mkdir -p $ROUT


##############################
### Colors
##############################

R CMD BATCH --no-save --no-restore "--args rwd='$RWD' out_dir='$RWD'" $RCODE/colors.R $ROUT/colors.Rout


##############################################################################
# Run DRIMSeq
##############################################################################


workers=2

simulation='drosophila_node_nonull'
count_method='kallistoprefiltered5'


for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull'
do
  for count_method in 'kallistoprefiltered5'
  do
    echo "${simulation}_${count_method}"

    R CMD BATCH --no-save --no-restore "--args rwd='$RWD/${simulation}' out_dir='drimseq_1_3_3/${count_method}' workers=${workers} count_method='${count_method}' one_way=TRUE coef_mode='optim' disp_moderation='trended' out_prefix='drimseq_oneway_trended_'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run.Rout
    tail $ROUT/sim5_drimseq_run.Rout

    R CMD BATCH --no-save --no-restore "--args rwd='$RWD/${simulation}' out_dir='drimseq_1_3_3/${count_method}' workers=${workers} count_method='${count_method}' one_way=TRUE coef_mode='optim' disp_moderation='none' out_prefix='drimseq_oneway_none_'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run.Rout
    tail $ROUT/sim5_drimseq_run.Rout

    R CMD BATCH --no-save --no-restore "--args rwd='$RWD/${simulation}' out_dir='drimseq_1_3_3/${count_method}' workers=${workers} count_method='${count_method}' one_way=FALSE coef_mode='optim' disp_moderation='trended' out_prefix='drimseq_reg_optim_trended_'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run.Rout
    tail $ROUT/sim5_drimseq_run.Rout

    R CMD BATCH --no-save --no-restore "--args rwd='$RWD/${simulation}' out_dir='drimseq_1_3_3/${count_method}' workers=${workers} count_method='${count_method}' one_way=FALSE coef_mode='optim' disp_moderation='none' out_prefix='drimseq_reg_optim_none_'" $RCODE/sim5_drimseq_run.R $ROUT/sim5_drimseq_run.Rout
    tail $ROUT/sim5_drimseq_run.Rout

  done
done


##############################################################################
### Comparison
##############################################################################

simulation='drosophila_node_nonull'
count_method='kallistoprefiltered5'


for simulation in 'drosophila_node_nonull' 'hsapiens_node_nonull'
do
  for count_method in 'kallistoprefiltered5'
  do

    echo "${simulation}_${count_method}"

    R CMD BATCH --no-save --no-restore "--args rwd='$RWD/${simulation}' out_dir='drimseq_1_3_3_comparison/${count_method}' workers=${workers} count_method='${count_method}' colors_path='$RWD/colors.Rdata' drimseq_res_dir='drimseq_1_3_3/${count_method}'" $RCODE/sim5_drimseq_comparison.R $ROUT/sim5_drimseq_comparison_run.Rout
    tail $ROUT/sim5_drimseq_comparison_run.Rout

  done
done










#
