##!/bin/sh
#$ -S /bin/sh
#current environment variables are used on you SGE jobs
#$ -V
#Use current working directory
#$ -cwd
#Merge the standard out and standard error to one file
#$ -j y
#$ -o $JOB_ID-output.txt 
#$ -t 1:10	   
#$ -l h_vmem=8G
#$ -q psy0cf20           # This is the computer queue. For high RAM jobs use: psy0cf20. 


cd /volume/DP_FEF/Analysis/28-Mar-2018

LOOPID='/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/loop_file.txt'
loopid=($(awk "NR==$SGE_TASK_ID" $LOOPID))

cd /volume/DP_FEF/Analysis/06-Apr-2018/temp_RHO
/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/dp_RHO_avg_5/for_testing/dp_RHO_avg_5 ${loopid[0]}

