## This script writes a sh script to send R jobs

writeSH <- function(nameSH,logName,mem,cmdScript,numJobs=1,numParJobs=1,numThreads=1)
{
  sink(nameSH)
  cat("#!/bin/sh\n")
  cat("#$ -S /bin/sh\n")
  cat("#$ -cwd\n")
  cat(paste("#$ -N",logName,"\n"))
  cat(paste0("#$ -t 1:",numJobs,"\n"))
  cat(paste("#$ -tc",numParJobs,"\n"))
  cat(paste0("#$ -l h_vmem=",mem,"\n"))
  cat(paste0("#$ -l tmem=",mem,"\n"))
  cat("#$ -j y\n")
  cat("#$ -R y\n")
  cat(paste0("#$ -pe smp ",numThreads,"\n"))
  cat("\n")
  cat("date\n")
  cat(paste0(cmdScript,"\n"))
  cat("date\n")
  cat("\n")
  sink()
}  


## writeSH("test.sh","test")
