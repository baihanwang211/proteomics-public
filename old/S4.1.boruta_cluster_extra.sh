#!/bin/bash
#SBATCH -A ckb.prj
#SBATCH -J S4.1.boruta_cluster_extra
#SBATCH -o S4.1.boruta_cluster_extra.out 
#SBATCH -e S4.1.boruta_cluster_extra.err 
#SBATCH -p short
#SBATCH -c 10

echo "------------------------------------------------" 

echo "Slurm Job ID: $SLURM_JOB_ID" 

echo "Run on host: "`hostname` 

echo "Operating system: "`uname -s` 

echo "Username: "`whoami` 

echo "Started at: "`date` 

echo "------------------------------------------------" 

module load R/4.2.1-foss-2022a

cd /gpfs3/well/ckb/users/dma206/proteomics/scripts

Rscript S4.1.boruta_cluster_extra.R