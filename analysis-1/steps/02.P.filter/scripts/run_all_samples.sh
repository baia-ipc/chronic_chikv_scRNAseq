#!/bin/bash

THISSCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export knit2html=$(readlink -f $THISSCRIPTDIR/../../../../scripts/knit2html)

SAMPLES="037-A 037-6m 229-A 372-A 267-6m 246-A 246-6m 227-A"
SAMPLES="$SAMPLES 227-6m 041-A 041-6m 262-6m 262-A 229-6m 217-6m"
SAMPLES="$SAMPLES 217-A 219-6m 266-A 266-6m 267-A"

n_threads=32

if [ "$2" != "" ]; then PSTR="prjpath=$1 libpath=$2"; else PSTR=""; fi

function run_sample { local SAMPLE=$1
  echo "Running sample $SAMPLE"
  ln -s -f filter_so.Rmd filter_so.$SAMPLE.Rmd
  $knit2html filter_so.$SAMPLE.Rmd sample=$SAMPLE $PSTR
  rm -f filter_so.$SAMPLE.Rmd
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

