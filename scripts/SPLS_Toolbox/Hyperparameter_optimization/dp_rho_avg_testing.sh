#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -o /volume/DP_FEF/Analysis/17-Apr-2018/$JOB_ID-output.txt   #output directory
#$ -j y
#$ -N DP_RHO_avg      # Name of the job
#$ -t 1-100:1       
#$ -l h_vmem=1G
#$ -q psy0cf20           # This is the computer queue. For high RAM jobs use: psy0cf20. 

/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/Hyperparameter_optimization/dp_RHO_avg_k_AR_DP/for_testing/dp_RHO_avg_k_AR_DP $SGE_TASK_ID
