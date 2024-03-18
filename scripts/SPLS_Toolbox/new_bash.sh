#!/bin/sh
# Tell the SGE that this is an array job, with "tasks" to be numbered 1 to 10000
#$ -t 1-10
#$ -tc 40
# When a single command in the array job is sent to a compute node,
# its task number is stored in the variable SGE_TASK_ID,
# so we can use the value of that variable to get the results we want:
#$ -i /volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/loop_file.$SGE_TASK_ID
#$ -o /volume/DP_FEF/Analysis/16-April-2018/output.$SGE_TASK_ID

/volume/DP_FEF/ScrFun/ScriptsRepository/SPLS_Toolbox/Hyperparameter_optimization/dp_RHO_avg_k_AR_DP/for_testing/dp_RHO_avg_k_AR_DP
