#!/bin/bash
#
#$ -cwd
#$ -o /volume/DP_FEF/Jobst_Bindung/Analysis/$JOB_ID-output-main.txt #output directory
#$ -j y
#$ -N SPLS_serial        # Name of the job
#$ -S /bin/bash
#$ -M david.popovic@med.uni-muenchen.de
#$ -m ae                # This is if you want emails when jobs abort or exit
#$ -l mem_free=2G       # This is how much memory you're requesting
#$ -l h_rss=2G          # This is the absolute memory limit
#$ -l h_stack=256M      # This is for matlab
#$ -q mitnvp1           # This is the computer queue. For high RAM jobs use: psy0cf20. 
cd /volume/DP_FEF/Jobst_Bindung/Analysis
/opt/matlab/R2015a/bin/matlab -nodisplay -nojvm -nosplash -r "run('/volume/DP_FEF/Jobst_Bindung/ScrFun/dp_aj_attachment.m')"