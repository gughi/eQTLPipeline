### this script plots the convergence of the PEER


setwd("/home/guelfi/eQTLPipeline/")

tmp <- read.delim(pipe("grep 'Delta bound:' testPEER/RNDMPEER18.log | cut -f6 -d' '"))
tmp <- gsub(",","",tmp[,1])
tmp <- as.numeric(tmp)

jpeg(paste0("plots/convergence.jpeg"), type="cairo")
plot(tmp[-c(1:50)], type="l", col="blue")
dev.off()

