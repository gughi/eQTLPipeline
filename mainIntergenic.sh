#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -N mainIntergenic
#$ -t 1:1
#$ -tc 1
#$ -j y
#$ -pe threaded 4

date
/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=/home/guelfi/eQTLPipeline/mainIntergenic.R
date
