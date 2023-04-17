#!/bin/bash
#SBATCH --job-name=VDJrecord
#SBATCH --time=24:00:00
#SBATCH --partition=rra
#SBATCH --qos=rra
#Notes below for how to estimate time based on what processes need to be run

#SBATCH --ntasks-per-node=8
#SBATCH --mem=40096

#SBATCH --output=/work/pi_gblanck/Arpan/NBL_WXS/SLURM_OUTPUT/out.t4VDJ.%j
#SBATCH --error=/work/pi_gblanck/Arpan/NBL_WXS/SLURM_OUTPUT/error.t4VDJ.%j

#SBATCH --mail-type=ALL
#SBATCH --mail-user=placeholder@aot.com

# MAKE SURE TO COPY sample.tsv TO NBL_BAMS_Results/final_csv

module add apps/python/3.8.5
python3 -m pip install localcider
python3 -m pip install lifelines
# -u makes it print stuff immediately instead of everything once it finishes.
python3 -u t4_run_VDJrecord_kon_dlbcl.py /work/pi_gblanck/Arpan/NBL_WXS/NBL_BAMS_Results/final_csv /work/pi_gblanck/Arpan/NBL_WXS/NBL_BAMS_Results/final_csv/sample.tsv /work/pi_gblanck/Arpan/NBL_WXS/scripts/threading_db
