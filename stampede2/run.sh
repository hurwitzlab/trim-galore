#!/bin/bash

module load tacc-singularity
module load launcher

PARAMRUN="$TACC_LAUNCHER_DIR/paramrun"

export LAUNCHER_PLUGIN_DIR="$TACC_LAUNCHER_DIR/plugins"
export LAUNCHER_WORKDIR="$PWD"
export LAUNCHER_RMI="SLURM"
export LAUNCHER_SCHED="interleaved"

set -u

IMG="/work/projects/singularity/TACC/biocontainers/trim-galore_0.4.5--pl5.22.0_0.img"
OUT_DIR="$PWD/trim-galore-out"

[[ ! -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"

if [[ ! -f "$IMG" ]]; then
    echo "Cannot find Singularity image \"$IMG\""
    exit 1
fi

QUERY=""
ARGS="-o $OUT_DIR"
while (($#)); do
    if [[ $1 == '-x' ]]; then
        QUERY="$QUERY $2"
        shift 2
    else
        ARGS="$ARGS $1"
        shift
    fi
done

if [[ -z "$QUERY" ]]; then
    echo "No -x QUERY"
    exit 1
fi

FILES_LIST=$(mktemp)
for QRY in $QUERY; do
    if [[ -f "$QRY" ]]; then
        echo "$QRY" >> "$FILES_LIST"
    elif [[ -d "$QRY" ]]; then
        find "$QRY" -type f >> "$FILES_LIST"
    else
        echo "QUERY \"$QRY\" neither file nor directory"
    fi
done

cat -n $FILES_LIST

NUM_FILES=$(wc -l "$FILES_LIST" | awk '{print $1}')

if [[ $NUM_FILES -lt 1 ]]; then
    echo "No good input files in \"$QUERY\""
    exit 1
fi

PARAM="$$.param"

i=0
while read -r FILE; do
    let i++
    printf "%3d: %s\n" $i $(basename "$FILE")
    echo "singularity exec $IMG trim_galore $ARGS $FILE" >> "$PARAM"
done < "$FILES_LIST"

export LAUNCHER_JOB_FILE="$PARAM"
$PARAMRUN

rm "$PARAM"

echo "Done, processed NUM_FILES \"$NUM_FILES\""
