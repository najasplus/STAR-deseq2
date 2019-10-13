#!/bin/bash
#$ -N STAR_GRCz11
#$ -pe parallel 24
#$ -l h_vmem=24G
#$ -o STAR_GRCz11_o.log
#$ -e STAR_GRCz11_e.log

./star_mapping.sh -s sample_description.txt -d '.' -o '--outFilterScoreMinOverLread 0.3'

