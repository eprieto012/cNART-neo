#!/bin/bash
#SBATCH -p highmem
#SBATCH --job-name=mixcr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=endikaprieto@vhio.net

set -euo pipefail

JAVA_HEAP="100g"

docker run --rm \
  --cpus "${SLURM_CPUS_PER_TASK}" \
  --memory "150g" \
  -e MI_LICENSE="E-VRKEGFHFHNBQIFCUJZOTZPNQISODMMGKRROTCOLBAOTTAOKA" \
  -e JAVA_TOOL_OPTIONS="-Xmx${JAVA_HEAP}" \
  -v /mnt/petasan_immuno/raw_data/itag/scSeq/Novogene/10x_Novogene_7/X208SC26017969-Z01-F001_01/01.RawData/ITAG_AGG_VDJ_01/:/raw:ro \
  -v /mnt/bioinfnas/immuno/eprieto/mixcr/:/work \
  ghcr.io/milaboratory/mixcr/mixcr:latest-develop \
  mixcr analyze 10x-sc-xcr-vdj-gemx-v3 \
    --species hsa \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --sample-sheet /work/TIL15_ITAG_AGG_S01/10x-sample-sheet.tsv \
    /raw/ITAG_AGG_VDJ_01-SCI7T002-SCI5T002_23GGJNLT3_S3_L006_R1_001.fastq.gz \
    /raw/ITAG_AGG_VDJ_01-SCI7T002-SCI5T002_23GGJNLT3_S3_L006_R2_001.fastq.gz \
    /work/TIL15_ITAG_AGG_S01/sample_result