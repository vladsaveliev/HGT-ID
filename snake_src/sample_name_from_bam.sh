#!/usr/bin/env bash

input=$1

samtools view -H ${input} | awk '{ \
    if ($1~/^@RG/) { \
        for (i; i<=NF; i++ ) { \
            if ($i ~ /^SM:/) { \
                col=i \
            } \
        }; \
        sub("SM:", "", $col); \
        print $col \
    } \
}' | head -1
