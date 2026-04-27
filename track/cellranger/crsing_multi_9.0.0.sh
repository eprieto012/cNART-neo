#!/bin/bash

#SBATCH -p highmem
#SBATCH --job-name=ITAG_AGG_S01	          # Job name
#SBATCH --mail-type=END,FAIL          	     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=endikaprieto@vhio.net    # Where to send mail
#SBATCH --ntasks=1  				     # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --mem=150gb                    	     # Job memory request
#SBATCH --time=24:00:00               		# Time limit hrs:min:sec
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --no-requeue

## command variables

id="ITAG_AGG_S01"                      #also docker folder name
projdir="/mnt/bioinfnas/immuno/eprieto/single_cell/ITAG_AGG_S01/"

singdir="/agrosmulti/"
referencedir="/mnt/petasan_general_bioinformatics_R/refs/sc10xgen/"
refsing="/refsing/"
transcriptome="refdata-gex-GRCh38-2024-A"
vdjref="refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
wdir="/currentdir/"
config="multiconfig.csv"

## main command

cmd="singularity run -H $PWD:$wdir -B $referencedir:$refsing -B $projdir/:$singdir -B /mnt/petasan_immuno/raw_data/itag/scSeq/Novogene/10x_Novogene_7:/rawdata cellranger9.0.0.sif \
     cellranger multi --id=$id \
     --csv=$wdir$config"
$cmd