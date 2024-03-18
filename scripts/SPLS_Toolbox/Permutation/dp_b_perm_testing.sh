#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -o /volume/DP_FEF/Analysis/17-Apr-2018/$JOB_ID-output.txt   #output directory
#$ -j y
#$ -N DP_PERM       # Name of the job
#$ -t 1-100:1    
#$ -tc 10
#$ -l h_vmem=1G
#$ -q psy0cf20           # This is the computer queue. For high RAM jobs use: psy0cf20. 

/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/Permutation/dp_perm_testing/for_testing/dp_perm_testing $SGE_TASK_ID

