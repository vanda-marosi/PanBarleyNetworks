#!/bin/bash
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=deseq2_panbarley
#SBATCH --output=/home/pgsb/vanda.marosi/logs/%x.%j.%a.out
#SBATCH --error=/home/pgsb/vanda.marosi/logs/%x.%j.%a.err
#SBATCH --array=1-20
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vanda.marosi@helmholtz-munich.de

#set -x

source ~/.bash_profile
source ~/.bashrc

# lookuptable for accession names
PROJECT_TABLE=/lustre/groups/pgsb/projects/panbarley_networks/00_geneTPM_peraccession/PanBaRT20_project_table.txt

name=$(cat $PROJECT_TABLE | head -n${SLURM_ARRAY_TASK_ID} | tail -n1)

OUTDIR=/lustre/groups/pgsb/projects/panbarley_networks/01_vsdTPMs_peraccession
CNTDIR=/lustre/groups/pgsb/projects/panbarley_networks/00_geneTPM_peraccession
METADIR=/lustre/groups/pgsb/projects/panbarley_networks/00_meta_peraccession

cnt_in="$CNTDIR/PanBaRT20_geneTPM_ort_filt_${name}.csv"
meta_in="$METADIR/${name}_meta.csv"
out="$OUTDIR/${name}_vsdTPM.csv"

date
conda deactivate
Rscript --vanilla /lustre/groups/pgsb/projects/panbarley_networks/scripts/02_varstabDESeq2.R $cnt_in $meta_in $out

echo "DESeq2 has run successfully."

# sbatch -N 1 /lustre/groups/pgsb/projects/panbarley_networks/scripts/run_02_deseq2_slurm.sh
