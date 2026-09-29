#!/bin/bash
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# Native GingerALE 3.0.2 ALE of published hippocampal-seed resting-state FC foci (healthy adults).
# Cluster-level FWE p<.05, cluster-forming p<.001, 1000 permutations, non-additive (Turkeltaub 2012), MNI_wb mask.
cd "$ROOT"
J="${GINGERALE_JAR:-/Applications/GingerALE.app/Contents/Java/GingerALE.jar}"
# GingerALE whole-brain MNI mask (extracted from the GingerALE jar; not redistributed here)
[ -f "$ROOT/04_results/gingerale/MNI_wb.nii" ] || (cd "$ROOT/04_results/gingerale" && unzip -o -q -j "$J" org/brainmap/image/MNI_wb.nii)
OUT=04_results/rsfc_ale/gingerale; mkdir -p $OUT
for f in "$@"; do
  b=$(basename $f .txt); cp $f $OUT/$b.txt
  (cd $OUT && java -Xmx8g -Djava.awt.headless=true -cp $J org.brainmap.meta.getALE2 $b.txt -mask=../../gingerale/MNI_wb.nii -nonadd -p=0.001 -perm=1000 -clust=0.05) > $OUT/$b.log 2>&1
  echo "$b done $(date)"
done
