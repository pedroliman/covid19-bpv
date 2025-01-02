#!/bin/bash
# Set Partition
#SBATCH --partition=short
# Request nodes
#SBATCH --nodes=1
# Request some processor cores
#SBATCH --ntasks=1
# Request memory 
#SBATCH --mem=16G
# Name of this job
#SBATCH --job-name=BVPCalcvar2000
# Output of this job, stderr and stdout are joined by default
# %x=job-name %j=jobid
#SBATCH --output=%x_%j.out
# Notify me via email
#SBATCH --mail-user=sarah.nowak@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array [1001-2000] ### Array index | %X number of simultaneous jobs
# 
# change to the directory where you submitted this script
cd /users/s/n/snowak/COVID-ABM
#
# your job execution follows:
echo "Starting sbatch script myscript.sh at:`date`"
# echo some slurm variables for fun
echo "  running host:    ${SLURMD_NODENAME}"
echo "  assigned nodes:  ${SLURM_JOB_NODELIST}"
echo "  partition used:  ${SLURM_JOB_PARTITION}"
echo "  jobid:           ${SLURM_JOBID}"
#
# load R
spack load r@3.6.3
spack load r-ggplot2@3.2.0
spack load r-dbplyr@1.4.2
spack load r-tidyr@0.8.3
spack load r-readr@1.3.1
spack load r-purrr@0.3.4
spack load r-tibble@2.1.3
spack load r-stringr@1.4.0
spack load r-forcats@0.4.0




Rscript BVP-dimless-parallel.R ${SLURM_ARRAY_TASK_ID} 






