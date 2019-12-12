IMG="/work/05066/imicrobe/singularity/trim-galore-0.6.5.img"

if [[ ! -e "$IMG" ]]; then
    echo "Missing Singularity image \"$IMG\""
    exit 1
fi

singularity exec $IMG run_trim_galore -o "trim-galore-out" ${QUERY} ${CONSIDER_ALREADY_TRIMMED} ${DONT_GZIP} ${GZIP} ${ILLUMINA} ${KEEP} ${LENGTH_1} ${LENGTH_2} ${NEXTERA} ${NON_DIRECTIONAL} ${PHRED64} ${RBBS} ${RETAIN_UNPAIRED} ${SMALL_RNA} ${TRIM} ${TRIM_N} ${ADAPTER} ${ADAPTER2} ${BASENAME} ${CLIP_R1} ${CLIP_R2} ${ERROR_RATE} ${HARDTRIM3} ${HARDTRIM5} ${LENGTH} ${MAX_LENGTH} ${MAX_N} ${NEXTSEQ} ${NUM_CONCURRENT_JOBS} ${NUM_HALT} ${QUALITY} ${STRINGENCY} ${THREE_PRIME_CLIP_R1} ${THREE_PRIME_CLIP_R2}

echo "Comments to kyclark@email.arizona.edu"
