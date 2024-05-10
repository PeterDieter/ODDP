#!/bin/bash

## SLURM Settings ################################################################
#SBATCH -t 0:10:00						# Wall clock time limit hrs:min:sec;
#															# 2 d = 48:40:00, 172800 ||
#															# 1 w = 168:40:00, 604800 ||
#															# 6:10 Std = 22000 ||
#															# 12:10 Std = 44000 ||
#															# 5:20 Std = 18000
#SBATCH -N 1									# (--nodes=<number>) Node count required for the job
#SBATCH -n 1									# (--ntasks=<number>) Number of tasks to be launched
#SBATCH --cpus-per-task=1			# CPU-Cores used per task
#SBATCH --ntasks-per-node=1		# Number of Tasks for the Job
#SBATCH --mem 4g							# Memory
##SBATCH --mem-per-cpu MEM		# Memory per allocated cpu core
#SBATCH -J dafd								# (--job-name=<jobname>) Job name
#SBATCH -A hpc-prf-pose				# (--account=<accountname>) Group
#SBATCH -p normal 						# (--partition <type>) Partition/queue in which to run job
##SBATCH --exclusive					# use compute node exclusively (not necessary)
#SBATCH --output log_dafd.out	# Standard output and error log
#SBATCH --error log_dafd.err	# Standard output and error log
#SBATCH --mail-type FAIL			# Mail events (NONE, BEGIN, END, FAIL, REQUEUE, ALL)
#SBATCH --mail-user psp@mail.upb.de
##SBATCH -q nocont						# priority in the queue ("express" uses your own user quota)

## silence core-allocation output in error-file
export SLURM_CPU_BIND=cores,quiet
export OMPI_MCA_hwloc_base_report_bindings=false

## clean module environment first
## reset module environment to system default (no output)
module reset > /dev/null 2>&1
## "Decide" which version by commenting out the desired version.
module load lib/gurobi/1101
## c++ compiler
module load compiler/GCC/13.2.0
## python version with venv
#module load lang/Python/3.9.6-GCCcore-11.2.0
#source ./venv/bin/activate

srun ./onlineAssignment --instance="instances/${1}.txt" --maxWaiting=${2} --AMethod=${3} --RMethod=${4} --${5} --a=${6} --beta=${7}

## deactivate venv
#deactivate
