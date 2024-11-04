#' @export
#'
print.jointLTSM <- function(x,...)
{
  
  if(class(x)!="jointLTSM") stop("Argument object should be a model from class 'jointLTSM'")
  
  cat("(Joint) Latent Time Shift Model, \n")
  cat("     fitted by maximum likelihood method", "\n")
  
  cl <- x$call
  cl$B <- NULL
  if(is.data.frame(cl$data))
  {
    cl$data <- NULL
    x$call$data <- NULL    
  }
  cat(" \n")
  dput(cl)
  cat(" \n")
  
  posfix <- eval(cl$posfix)
  
  cat("Statistical Model:", "\n")
  cat(paste("     Dataset:", as.character(as.expression(x$call$data))),"\n")
  cat(paste("     Number of subjects:", x$ns),"\n")
  
  cat("     Number of observations:", paste(x$N[[2]]),"\n")
  cat(paste("     Number of parameters:", length(x$best))," \n")
  if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")
  
  
  nbevt <- x$N[[3]]
  if(nbevt>0)
  {
    nprisq <- rep(NA,nbevt)
    nrisq <- rep(NA,nbevt)
    for(ke in 1:nbevt)
    {
      if(x$typrisq[ke]==1) nprisq[ke] <- x$nz[ke]-1
      if(x$typrisq[ke]==2) nprisq[ke] <- 2
      if(x$typrisq[ke]==3) nprisq[ke] <- x$nz[ke]+2
      
      
      cat(paste("     Event",ke,": \n"))
      cat(paste("        Number of events: ", x$nevent[ke],"\n",sep=""))
      
      if (x$typrisq[ke]==2)
      {
        cat("        Weibull baseline risk function \n")
      }
      if (x$typrisq[ke]==1)
      {
        cat("        Piecewise constant baseline risk function with nodes \n")
        cat("        ",x$hazardnodes[1:x$nz[ke],ke]," \n")
      }
      if (x$typrisq[ke]==3)
      {
        cat("        M-splines constant baseline risk function with nodes \n")
        cat("        ",x$hazardnodes[1:x$nz[ke],ke]," \n")
      }
      
    }
    cat("        Association: \n")
    cat("           Time-shift \n")
    for(yk in 1:x$N[[1]])
      cat(paste("           Outcome ",x$Names$Ynames[yk], " : ", ifelse(x$sharedtype[ke,yk]==0,"shared random effects","current value"),"\n",sep=""))
    
  }
  
  
  ny <- x$N[[1]]
  Ynames <- x$Names$Yname
  cat("     Link functions:  \n")
  for (yk in 1:ny)
  {
    if (x$linktype[yk]==0)
    {
      cat("        Identity for",Ynames[yk]," \n")
    }
    if (x$linktype[yk]==1)
    {
      cat("        Linear for",Ynames[yk], "\n")
    }
    if (x$linktype[yk]==2)
    {
      cat("        Quadratic I-splines with nodes", x$linknodes[1:x$nbzitr[yk],yk]," for ",Ynames[yk], "\n")
    }
  }
  
  
  
  cat(" \n")
  cat("Iteration process:", "\n")
  
  if(x$conv==1) cat("     Convergence criteria satisfied")
  if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
  if(x$conv==3) cat("     Convergence with restrained Hessian matrix")
  if(x$conv==4|x$conv==12)
  {
    cat("     The program stopped abnormally. No results can be displayed.\n")
  }
  else
  {
    cat(" \n")
    cat("     Number of iterations: ", x$niter, "\n")
    cat("     Convergence criteria: parameters=", signif(x$gconv[1],2), "\n")
    cat("                         : likelihood=", signif(x$gconv[2],2), "\n")
    cat("                         : second derivatives=", signif(x$gconv[3],2), "\n")
    cat(" \n")
    cat("Goodness-of-fit statistics:", "\n")
    cat(paste("     maximum log-likelihood:", round(x$loglik,2))," \n")
    cat(paste("     AIC:", round(x$AIC,2))," \n")
    cat(paste("     BIC:", round(x$BIC,2))," \n")
    cat(" \n")
    
  }
}
