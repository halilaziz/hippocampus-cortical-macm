#!/bin/bash
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# Compute ALE maps natively in GingerALE 3.0.2 (non-additive) for each pooled/neurosynth seed, for voxelwise
# equivalence testing against NiMARE. No permutation thresholding (too slow for >1000 experiments in GingerALE).
cd "$ROOT/04_results/gingerale"
J="${GINGERALE_JAR:-/Applications/GingerALE.app/Contents/Java/GingerALE.jar}"
# GingerALE whole-brain MNI mask (extracted from the GingerALE jar; not redistributed here)
[ -f "$ROOT/04_results/gingerale/MNI_wb.nii" ] || (cd "$ROOT/04_results/gingerale" && unzip -o -q -j "$J" org/brainmap/image/MNI_wb.nii)
mkdir -p out
for f in foci/neurosynth_*HIPP.txt foci/pooled_*HIPP.txt; do
  b=$(basename $f .txt)
  [ -f out/${b}_ALE.nii ] && continue
  cp $f out/$b.txt
  (cd out && java -Xmx16g -Djava.awt.headless=true -cp $J org.brainmap.meta.getALE2 $b.txt -mask=../MNI_wb.nii -nonadd -noPVal) > out/$b.log 2>&1
  echo "$b done $(date)"
done
