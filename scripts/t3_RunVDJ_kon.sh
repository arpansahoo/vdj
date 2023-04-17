#!/bin/bash

#SBATCH --job-name=VDJ_Job
#SBATCH --time=72:00:00
#SBATCH --output=/work/pi_gblanck/Arpan/NBL_WXS/SLURM_OUTPUT/out.RunVDJ.%A_%a
#SBATCH --error=/work/pi_gblanck/Arpan/NBL_WXS/SLURM_OUTPUT/err.RunVDJ.%A_%a
#SBATCH --array=0-20
#SBATCH --mem=100250
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1

#SBATCH --mail-type=ALL
#SBATCH --mail-user=placeholder@aot.com

# need to change array=0-numjobs depending on console output from t3_set_task_items.py

module load apps/python/3.8.5
IFS='#' read -a receptors_to_do < /work/pi_gblanck/Arpan/NBL_WXS/scripts/t3_rectodo.txt
#receptors_to_do=(TRA TRB TRD TRG IGH IGK IGL TRA_UM TRB_UM TRD_UM TRG_UM IGH_UM IGK_UM IGL_UM)

echo "array task id:" $SLURM_ARRAY_TASK_ID
echo "receptor:" ${receptors_to_do[$SLURM_ARRAY_TASK_ID]}
python3 -u "/work/pi_gblanck/Arpan/NBL_WXS/scripts/t3_findvdjum.py" ${receptors_to_do[$SLURM_ARRAY_TASK_ID]} "/work/pi_gblanck/Arpan/NBL_WXS/NBL_BAMS_Results/" "/work/pi_gblanck/Arpan/NBL_WXS/scripts/t3_vdjdb/"

echo 'done'
