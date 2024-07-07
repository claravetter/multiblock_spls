#!/bin/sh
# Tell the SGE that this is an array job, with "tasks" to be numbered 1 to 10000
#$ -t 1-10000
# When a single command in the array job is sent to a compute node,
# its task number is stored in the variable SGE_TASK_ID,
# so we can use the value of that variable to get the results we want:
~/programs/program -i ~/data/input.$SGE_TASK_ID -o ~/results/output.$SGE_TASK_ID