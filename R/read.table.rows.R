read.table.rows <- 
function(filename, keepRows, sep="\t", header=T,rnames=1, ...)
{
  
  tmpf <- tempfile()
  write.table(data.frame(sort(keepRows)), file=tmpf, row.names=F, col.names=F, quote=F )
  
  out  <- read.table( pipe( paste("fgrep -w -f", tmpf, filename) ), sep=sep, header=F, row.names=rnames, ...)
  notFound <- setdiff( keepRows, rownames(out) )
  if( length(notFound) > 0) cat( "Warning:", length(notFound), "requested rows not found.\n")
  
  if(header) colnames(out) <- scan(filename, sep=sep, what=character(0), nlines=1, quiet=T)[-1]
  
  unlink(tmpf); rm(tmpf); gc()
  return(out)
}
# Function from Adai
