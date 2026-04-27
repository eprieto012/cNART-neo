#!/bin/bash

#SBATCH -p highmem
#SBATCH --job-name=mixcr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=endikaprieto@vhio.net
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=100gb
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

docker run --rm \
  -e MI_LICENSE="E-VRKEGFHFHNBQIFCUJZOTZPNQISODMMGKRROTCOLBAOTTAOKA" \
  -v /mnt/petasan_immuno/raw_data/itag/scSeq/Novogene/10x_Novogene_7/X208SC26017969-Z01-F001_01/01.RawData/ITAG_AGG_VDJ_02/:/raw:ro \
  -v /mnt/bioinfnas/immuno/eprieto/mixcr/:/work \
  ghcr.io/milaboratory/mixcr/mixcr:latest-develop \
  mixcr analyze 10x-sc-xcr-vdj-gemx-v3 \
    --species hsa \
    --sample-sheet /work/TIL15/10x-sample-sheet.tsv \
    /raw/ITAG_AGG_VDJ_02-SCI7T014-SCI5T014_23GGJNLT3_S4_L006_R1_001.fastq.gz \
    /raw/ITAG_AGG_VDJ_02-SCI7T014-SCI5T014_23GGJNLT3_S4_L006_R2_001.fastq.gz \
    /work/TIL15/sample_result