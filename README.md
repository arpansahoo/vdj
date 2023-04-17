# VDJ mining software 
This software, created in the Blanck Lab at the USF Morsani College of Medicine, recovers TCR V(D)J recombination reads from sequencing data and identifies CDR3 regions within those reads. We also include code for physico-chemical analysis of the CDR3s. The latest version of the software also includes a `threading` method in `t4_vdjrecord.py` to stitch V- and J- segments to the CDR3s, creating full-length V-CDR3-J sequences that represent the entire TCR variable region.

## Instructions
- Run `sbatch gdcslice.sh`
  - This will download bam slices from your manifest.
  - Get a manifest (for the WXS or RNA-seq files you're interested in) and token file beforehand from GDC
- Run `sbatch t2_Module_Search_IgTcR_header.sh`
  - First step in processing the bams
  - Edit the file paths in this script at the top and bottom
  - May need to edit the file paths in `t2_Module_Search_IgTcRFix.sh`
  - Make sure GNU parallel is installed
- Run `python t3_set_task_items.py`
  - Edit the path at the top of the file to match the results folder generated from the previous step
  - Note the array setting printed out to put into the next slurm config file, `t3_Run_VDJ.sh`
- Run `sbatch t3_RunVDJ.sh`
  - The t3 module has multiple tasks, but most importantly, it finds the matching V/D/J/CDR3 sequence for each read 
  - Edit the file paths in this script
  - Remember to place the vdjdb folder into the directory that you indicate
- Run `sbatch t4_pre.sh`
  - Edit filepaths in `t4_pre.sh` and `t4_pre.py`
- Run `sbatch t4_start_VDJRecord.sh`
  - Edit cancer in `t4_run_VDJrecord.py`, edit filepaths in `t4_start_VDJrecord.sh`
  - Before running t4, download `sample.tsv` from GDC, and put it in the results from t4_pre
  
