#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -N compress
#$ -t 1:11
#$ -tc 1
#$ -j y
#$ -pe threaded 1
#$ -q highmem.q

date
/home/guelfi/softwares/R-3.2.0/bin/R --vanilla --file=/home/guelfi/eQTLPipeline/mainIntergenic.R
date
