#!/bin/bash
#
#$ -cwd
#$ -o /volume/DP_FEF/Analysis/06-Apr-2018/temp_RHO/ #output directory
#$ -j y
#$ -N RHO_avg_1_10        # Name of the job
#$ -S /bin/bash
#$ -M david.popovic@med.uni-muenchen.de
#$ -m ae                # This is if you want emails when jobs abort or exit
#$ -l mem_free=2G       # This is how much memory you're requesting
#$ -l h_rss=2G          # This is the absolute memory limit
#$ -l h_stack=256M      # This is for matlab
#$ -t 1:11	   
#$ -q psy0cf20          # This is the computer queue. For high RAM jobs use: psy0cf20. 


LOOPID='/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/loop_file.txt'
loopid=($(awk "NR==$SGE_TASK_ID" $LOOPID))

cd /volume/DP_FEF/Analysis/06-Apr-2018/temp_RHO/
/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/dp_RHO_avg_6/for_testing/dp_RHO_avg_6 ${loopid[0]}

