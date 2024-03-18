##!/bin/bash
#$ -S /bin/sh
#current environment variables are used on you SGE jobs
#$ -V
#Use current working directory
#$ -cwd
#Merge the standard out and standard error to one file
#$ -j y
#$ -o $JOB_ID-output.txt 
#$ -t 1:1	   
#$ -l h_vmem=8G
#$ -q psy0cf20           # This is the computer queue. For high RAM jobs use: psy0cf20. 


LOOPID='/volume/DP_FEF/Analysis/16-April-2018/missing_ID.txt'
loopid=($(awk "NR==$SGE_TASK_ID" $LOOPID))

/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/Hyperparameter_optimization/dp_RHO_avg_k_AR_DP/for_testing/dp_RHO_avg_k_AR_DP ${loopid[0]}

