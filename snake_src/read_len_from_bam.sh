#!/usr/bin/env bash

input=$1

samtools view $input | \
head -n10000 | \
awk '{ print length($10) }' | \
awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'
