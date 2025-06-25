#!/bin/bash

THISSCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export knit2html=$(readlink -f $THISSCRIPTDIR/../../../../scripts/knit2html)

SAMPLES="037 S3-A3 S15-A7 S16-A8 271 S5-A5 S13-B3"
SAMPLES="$SAMPLES S14-B4 229 S4-A4 S1-A1 S2-A2 S7-A9 S8-A10 S9-A11"
SAMPLES="$SAMPLES S10-A12 S11-B1 S12-B2 372 S6-A6"

n_threads=32

if [ "$2" != "" ]; then PSTR="prjpath=$1 libpath=$2"; else PSTR=""; fi

function run_sample { local SAMPLE=$1
  echo "Running sample $SAMPLE"
  ln -s -f create_so.Rmd create_so.$SAMPLE.Rmd
  $knit2html create_so.$SAMPLE.Rmd sample=$SAMPLE $PSTR
  rm -f create_so.$SAMPLE.Rmd
  echo "Finished sample $SAMPLE"
}
export -f run_sample

NSAMPLES=$(echo $SAMPLES | wc -w)
if command -v parallel >/dev/null 2>&1; then
  echo "Running $NSAMPLES samples in parallel using $n_threads threads"
  parallel --line-buffer -j $n_threads run_sample ::: $SAMPLES
else
  echo "Running $NSAMPLES samples in serial as GNU parallel is not installed"
  for SAMPLE in $SAMPLES; do run_sample $SAMPLE; done
fi

