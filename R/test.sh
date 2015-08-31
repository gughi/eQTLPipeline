#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -N test 
#$ -t 1: 1 
#$ -tc 1#$ -j y
#$ -pe threaded 1 
