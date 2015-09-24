qq.plot <- function( x, alpha=0.05, datatype="pvalue", scaletype="pvalue", df=1, plot.concentration.band=TRUE, one.sided=FALSE, frac=1, iplot=NULL, print=FALSE, xat=NULL, yat=NULL, main=NULL, xlab=NULL, ylab=NULL, pch="x", cex=0.5, col="black", ... )
{
  pname <- paste(deparse(substitute(x), 500), collapse="\n")   #Name of vector passed as "x" to be used in plot title etc.
  if (!is.null(iplot)) frac=1            #Forces frac=1 if "iplot" used to chose points to plot
  #Some validity checks on x
  if (!is.numeric(x))
    stop("'x' must be numeric")
  nmissing = sum(is.na(x))
  x <- x[ !is.na(x) ]                 #To deal with missing data values (these don't get plotted)
  if ((datatype=="pvalue")&((min(x)<0)|(max(x)>1)))
    stop("'x' must be >=0 and <=1")
  if ((datatype=="chisq")&(min(x)<0))
    stop("'x' must be >=0")
  nzero = sum(x==0)
  #Some warnings on missing values (and, if pvalues, on x=0) 
  if (nmissing>0)
    warning(nmissing, " missing values (excluded)")
  if ((nzero>0)&(datatype=="pvalue")) {
    warning(nzero, " zero values (plotted with same value as lowest non-zero p-value)")
    x[x==0] <- min(x[x>0])
  }
  if (datatype=="stdnorm") {df=0; scaletype="ordinal"}
  n <- length(x)
  starti = floor((n-1)*(1-frac)) +1         			#i for the first sorted datapoint to be plotted.
  lena = n-starti+1			               			#Number of datapoints to be plotted
  if (!is.null(iplot)) a2=iplot else a2=(1:lena)          #indices to be plotted
  b <- n+1-a2                                             #indices used in determining concentration band
  #Find E and O under relevant inv. chisq transformation
  if ((df==2)&(datatype!="stdnorm")) {                                  #short-cut for df=2 (chisq or pval data): use -2log-transformed expected U(0,1) order statistics (high values first)
    E <- -2*log(a2/(n+1))                                               
    if (datatype=="pvalue") O <- -2*log(sort(x)[a2])                    #Note obs data no need to transform if already chisq or z value (high values first)
  } else {
    if (datatype=="stdnorm") E <- qnorm(a2/(n+1),lower.tail=FALSE)              #invnorm-transformed expected U(0,1) order statistics (put high scores first)
    if (datatype!="stdnorm") E <- qchisq(a2/(n+1),df=df,lower.tail=FALSE)       #invchisq-transformed expected U(0,1) order statistics (put high scores first)
    if (datatype=="pvalue") O <- qchisq(sort(x)[a2],df=df,lower.tail=FALSE)     #Take lowest pvalues, transform to chisq (highest/most interesting values first)
  }
  if (datatype!="pvalue") O <- sort(x, decreasing=TRUE)[a2]                     #Sort x (chisq or norm), highest (most interesting) values first
  #Derive "pretty" tick places for log10 p-value scale, if necessary
  #Note that by this stage, O/E will either contain chisq-scale or normal-scale values, and both are sorted
  if (scaletype=="pvalue") {                                              #Note scaletype forced="quantile" for stdnorm data, so here all data is on chisq scale
    if (!is.null(xat)) x4Lx=xat  else x4Lx = pretty( -log10( pchisq(c(E[1],E[length(E)]),df=df,lower.tail=FALSE) ) )
    if (!is.null(yat)) y4Ly=yat  else y4Ly = pretty( -log10( pchisq(c(O[1],O[length(O)]),df=df,lower.tail=FALSE) ) )
    xnums = qchisq(10^-x4Lx,df=df,lower.tail=FALSE)                       #Get same locations on actual chisq scale
    ynums = qchisq(10^-y4Ly,df=df,lower.tail=FALSE)                       #Get same locations on actual chisq scale
    Lx <- parse( text=paste("10^-",x4Lx,sep="") )
    Ly <- parse( text=paste("10^-",y4Ly,sep="") )
  } else {                                                                #"Else" covers both chisq and stdnorm-scaled data
    if (!is.null(xat)) xnums=xat else xnums=pretty(c(E[1],E[length(E)]))
    if (!is.null(yat)) ynums=yat else ynums=pretty(c(O[1],O[length(O)]))
    Lx <- parse( text=as.character(xnums) )
    Ly <- parse( text=as.character(ynums) )
  }
  #Do Q-Q plot
  if (is.null(main)) {
    if (datatype=="stdnorm") main=paste("Q-Q plot (on stdnorm) of " ,pname, sep="")
    else main=paste("Q-Q plot (on chisq[",df,"]) of " ,pname, sep="")
  }
  if (is.null(xlab)) {
    if (scaletype=="pvalue") xlab="Expected p-value"
    else xlab="Expected quantile"
  }
  if (is.null(ylab)) {
    if (scaletype=="pvalue") ylab="Observed p-value"
    else ylab="Observed quantile"
  }
  plot( c(E[1],E[length(E)]), c(O[1],O[length(O)]), main = main, xlab = xlab, ylab = ylab, type = "n", xaxt = "n", yaxt = "n", ... )       #Just plots the outside box
  axis(1, at=xnums, labels=Lx )
  axis(2, at=ynums, labels=Ly )
  if (plot.concentration.band==TRUE) {        #Note that conc band won't draw if x has too many datapoints
    if (one.sided==FALSE) {
      upper <- qbeta( 1-alpha/2, a2, b )      #Exp. upper CL for 'a'th U(0,1) order statistic (becomes 'lower')
      lower <- qbeta( alpha/2, a2, b )        #Exp. lower CL for 'a'th U(0,1) order statistic (becomes 'upper')
    } else {
      upper <- rep(1,length(E))                    #Exp. upper CL for 'a'th U(0,1) order statistic (becomes 'lower')
      lower <- qbeta( alpha, a2, b )          #Exp. lower CL for 'a'th U(0,1) order statistic (becomes 'upper')
    }
    if (df==2) {
      polygon( c( E, rev(E) ), c( -2*log(upper), rev(E) ), col="grey", border = NA )  #'lower' band after trans
      polygon( c( E, rev(E) ), c( -2*log(lower), rev(E) ), col="grey", border = NA )  #'upper' band after trans
    } else {
      if (datatype=="stdnorm") {
        polygon( c( E, rev(E) ), c( qnorm(upper,lower.tail=FALSE), rev(E) ), col="grey", border = NA )
        polygon( c( E, rev(E) ), c( qnorm(lower,lower.tail=FALSE), rev(E) ), col="grey", border = NA )
      } else {
        polygon( c( E, rev(E) ), c( qchisq(upper,df=df,lower.tail=FALSE), rev(E) ), col="grey", border = NA ) #'lower' band
        polygon( c( E, rev(E) ), c( qchisq(lower,df=df,lower.tail=FALSE), rev(E) ), col="grey", border = NA ) #'upper' band
      }
    }
  }
  abline( 0, 1, col="red" )                                          #plot 1:1 line
  abline(h=ynums, v=xnums, col="lightgray", lty="dotted")            #plot grid
  points( E, O, pch=pch, cex=cex, col=col )                          #Finally, plot points
  if (print==TRUE) return( data.frame( O=O, E=E ) )
}

