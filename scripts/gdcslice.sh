#!/bin/bash
#SBATCH --job-name=GDC-SliceDownload
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mem-per-cpu=2048
#SBATCH --output=/work/pi_gblanck/Arpan/NBL_WXS/SLURM_OUTPUT/out.GDC_SLICE.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=placeholder@aot.com

# change this number to the total number of files in your manifest, ie the number 443; do not change the %100
#SBATCH --array=1-443%100

# change this to your base path
basepath="/work/pi_gblanck/Arpan/NBL_WXS"

# include manifest file name as argument for working with clusters
# comment out if not working with clusters
#manifest=$1

# path to manifest, change this to your manifest location
# for small manifest files
PathToManifest="${basepath}/manifest.txt"

# for clusters
#PathToManifest="${basepath}/${manifest}"
#echo $PathToManifest

#path to download token, change this to your token location
Token="${basepath}/token.txt"

# change this path to the path you want the BAM files to go to
OutputFolder="${basepath}/BAMS"

mapfile -t myArray < $PathToManifest

NumberOfBams=$((${#myArray[@]} - 1))
echo "There are" $NumberOfBams "bam files"
echo "$InputString"
InputString=${myArray[$SLURM_ARRAY_TASK_ID]}
ID=$(cut -d' ' -f1 <<< $InputString)
NAME=$(cut -d' ' -f2 <<< $InputString)

APItext="https://api.gdc.cancer.gov/slicing/view/$ID?region=chr14:21521904-22652132&region=chr7:142289011-142913287&region=chr7:38240024-39368055&region=chr14:105486437-106879844&region=chr2:88857361-90238368&region=chr22:22026076-22922913&region=chr6:29844528-33100696&region=unmapped"
#                                                             TRA / TRD                       TRB                              TRG                          IGH                                 IGK                               IGL                            HLA
token=$(<$Token)
curl --header "X-Auth-Token: $token" $APItext --output $OutputFolder/sliced_$NAME