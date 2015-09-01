## This script writes a sh script to send R jobs

writeSH <- function(nameSH,logName, script,numJobs=1,numParJobs=1,numThreads=1)
{
  sink(nameSH)
  cat("#!/bin/sh\n")
  cat("#$ -S /bin/sh\n")
  cat("#$ -cwd\n")
  cat(paste("#$ -N",logName,"\n"))
  cat(paste("#$ -t 1:",numJobs,"\n"))
  cat(paste("#$ -tc",numParJobs,"\n"))
  cat("#$ -j y\n")
  cat(paste("#$ -pe threaded",numThreads,"\n"))
  cat("\n")
  cat(paste(script,"\n"))
  cat("\n")
  sink()
}  
  
  
writeSH("test.sh","test")

