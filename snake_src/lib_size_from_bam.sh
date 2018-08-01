#!/usr/bin/env bash

input=$1

samtools view -q 20 -f2 ${input} | \
awk '$6 !~ /N/' | \
awk '$7 ~ /=/' | \
cut -f9 | awk '$1<1000' | \
head -10000 | \
awk '{if ($1<0) {$1=-$1} else {$1=$1} sum+=$1;} END { printf("%3.0f", sum/NR) }'

