#!/bin/bash
#SBATCH -A ckb.prj
#SBATCH -J boruta_dedup
#SBATCH -p short

echo "------------------------------------------------" 

echo "Slurm Job ID: $SLURM_JOB_ID" 

echo "Run on host: "`hostname` 

echo "Operating system: "`uname -s` 

echo "Username: "`whoami` 

echo "Started at: "`date` 

echo "------------------------------------------------" 

module load R/4.2.1-foss-2022a

cd /gpfs3/well/ckb/users/dma206/proteomics/scripts

Rscript correlation_explain_deduplicated_cluster.R