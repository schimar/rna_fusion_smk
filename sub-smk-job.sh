#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# sub-smk-job.sh — run the wts (rna_fusion) Snakemake workflow inside the
# central Docker image.
#
# Usage:
#   ./sub-smk-job.sh -n --config runid=/path/to/run bcldir=/mnt/illumina/runfolder analysis_path=analysis_rna
#   ./sub-smk-job.sh -j 50 --config runid=/path/to/run bcldir=/mnt/illumina/runfolder analysis_path=analysis_rna
#
# Any additional arguments are passed straight through to `snakemake`.
# Set RNA_FUSION_IMAGE to override the default image tag.
# ---------------------------------------------------------------------------
set -euo pipefail

IMG="${RNA_FUSION_IMAGE:-ukwgenommedizin/rna_fusion:1.0.1}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec docker run --rm \
    -u "$(id -u):$(id -g)" \
    -e XDG_CACHE_HOME=/tmp/smk_cache \
    -v /tmp:/tmp \
    -v "/home:/home" \
    -v "/mnt:/mnt" \
    -v "$REPO:$REPO" \
    -w "$REPO/workflow" \
    "$IMG" \
    snakemake "$@"