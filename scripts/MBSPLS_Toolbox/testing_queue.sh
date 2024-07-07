##!/bin/bash
#$ -S /bin/sh
#current environment variables are used on you SGE jobs
#$ -V
#Use current working directory
#$ -cwd
#Merge the standard out and standard error to one file
#$ -j y
#$ -o $JOB_ID-output.txt    
#$ -l h_vmem=8G
#$ -pe 5
#$ -q psy0cf20           # This is the computer queue. For high RAM jobs use: psy0cf20. 


LOOPID='/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/test.txt'
loopid=($(awk "NR==$SGE_TASK_ID" $LOOPID))

/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/Hyperparameter_optimization/dp_RHO_avg_k_AR_DP/for_testing/dp_RHO_avg_k_AR_DP ${loopid[0]}

