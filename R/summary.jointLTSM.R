#' Summary of a (joint) latent time shift model
#' 
#' This function provides a summary of model estimated with the \code{jointLTSM} function.
#' 
#' @param object an object inheriting from class \code{jointLTSM} for a joint latent time shift model.
#' @param ... further arguments to be passed to or from other methods. They are ignored in this function.
#' 
#' @return The function prints a summary of the model estimates.
#' 
#' @author Tiphaine Saulnier, Cecile Proust-Lima
#' 
#' @export
#'
summary.jointLTSM <- function(object,...)
{
    x <- object
    
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


            cat("Maximum Likelihood Estimates:", "\n")
            cat(" \n")
            
            nobs <- x$N[[2]]
            nbevt <- x$N[[3]]
            nrisqtot <- x$N[[4]]
            nvarxevt <- x$N[[5]]
            nasso <- x$N[[6]]
            nef.s <- x$N[[7]]
            nvc.s <- x$N[[8]]
            nef <- x$N[[9]]
            nea <- x$N[[10]]
            nvc <- x$N[[11]]
            ntr <- x$N[[12]]
            nvc.err <- x$N[[13]]
            idlink <- x$linktype
            idmodel <- x$modeltype
            NPM <- length(x$best)
            
            se <- rep(NA,NPM)
            if (x$conv==1 | x$conv==3)
                {
                    ##recuperation des indices de V
                    id <- 1:NPM
                    indice <- id*(id+1)/2
                    se <- sqrt(x$V[indice])
                    se[which(is.na(se))] <- 1
                    wald <- x$best/se
                    pwald <- 1-pchisq(wald**2,1)
                    coef <- x$best
                }
            else
                {
                    se <- NA
                    wald <- NA
                    pwald <- NA
                    coef <- x$best

                    sech <- rep(NA,length(coef))
                    waldch <- rep(NA,length(coef))
                    pwaldch <- rep(NA,length(coef))
                }

            # valeur absolue des std.err
            for(k in 1:ny)
              coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+(k-1)*(sum(nef[1:(k-1)])+sum(nvc[1:(k-1)])+sum(ntr[1:(k-1)])+(k-1))+nef[k]+nvc[k]+ntr[k]+1] <- abs(coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+(k-1)*(sum(nef[1:(k-1)])+sum(nvc[1:(k-1)])+sum(ntr[1:(k-1)])+(k-1))+nef[k]+nvc[k]+ntr[k]+1])
              
            
            if(x$conv!=2)
                {
                    coefch <- format(as.numeric(sprintf("%.5f",coef)),nsmall=5,scientific=FALSE)
                    sech <- format(as.numeric(sprintf("%.5f",se)),nsmall=5,scientific=FALSE)
                    waldch <- format(as.numeric(sprintf("%.3f",wald)),nsmall=3,scientific=FALSE)
                    pwaldch <- format(as.numeric(sprintf("%.5f",pwald)),nsmall=5,scientific=FALSE)
                }
            else
                {
                    coefch <- format(as.numeric(sprintf("%.5f",coef)),nsmall=5,scientific=FALSE)
                }
            
            
            if(length(posfix))
                {
                    coefch[posfix] <- paste(coefch[posfix],"*",sep="")
                    sech[posfix] <- ""
                    waldch[posfix] <- ""
                    pwaldch[posfix] <- ""
                }

            ## fct pr determiner la longueur max d'une chaine de caracteres
            ## (avec gestion des NA)
            maxchar <- function(x)
                {
                    xx <- na.omit(x)
                    if(length(xx))
                        {
                            res <- max(nchar(xx))
                        }
                    else
                        {
                            res <- 2
                        }
                    return(res)
                }

#browser()

            if(nbevt>0)
            {
                cat("Survival model:\n" )
                cat("\n")
                cat("Parameters in the proportional hazard model:\n" )
                cat("\n")

                tmp <- cbind(coefch[1:(nrisqtot+nvarxevt+nasso)],
                             sech[1:(nrisqtot+nvarxevt+nasso)],
                             waldch[1:(nrisqtot+nvarxevt+nasso)],
                             pwaldch[1:(nrisqtot+nvarxevt+nasso)])
                maxch <- apply(tmp,2,maxchar)
                if(any(1:(nrisqtot+nvarxevt+nasso) %in% posfix)) maxch[1] <- maxch[1]-1
                dimnames(tmp) <- list(names(coef)[1:(nrisqtot+nvarxevt+nasso)],
                                      c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                        paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                        paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                        paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
                cat("\n")
                print(tmp,quote=FALSE,na.print="")
                cat("\n")
            }

            ##############################
            
            cat("-------------------------------------------------------------------\n" )
            cat("Time shift:\n" )
            cat(" \n")
            
            tmp <- cbind(coefch[nrisqtot+nvarxevt+nasso+1:nef.s],
                         sech[nrisqtot+nvarxevt+nasso+1:nef.s],
                         waldch[nrisqtot+nvarxevt+nasso+1:nef.s],
                         pwaldch[nrisqtot+nvarxevt+nasso+1:nef.s])
            tmp <- rbind(tmp,
                         c(coefch[nrisqtot+nvarxevt+nasso+nef.s+1:nvc.s],NA,NA,NA))
            names(coef)[nrisqtot+nvarxevt+nasso+nef.s+1:nvc.s] <- "s. residual variance"
            
            maxch <- apply(tmp,2,maxchar)
            if(any(c(nrisqtot+nvarxevt+nasso+1:(nef.s+nvc.s)) %in% posfix)) maxch[1] <- maxch[1]-1
            dimnames(tmp) <- list(names(coef)[nrisqtot+nvarxevt+nasso+1:(nef.s+nvc.s)],
                                  c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                    paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                    paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                    paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
            cat("\n")
            print(tmp,quote=FALSE,na.print="")
            cat("\n")
            
            ##########
            
            cat("-------------------------------------------------------------------\n" )
            cat("Outcome SubModel(s):\n" )
            cat("\n")
            
            sumPrmK = 0
            for(k in 1:ny){
              
              cat("\n")
              if(k>1) cat("---------------------------------\n" )
              cat(paste("  Outcome ",k," : ",x$Names$Ynames[k],"\n",sep=""))
              cat("\n")
              
              
              # trajectory
              cat(paste("    trajectory: ",ifelse(idmodel[k]==0,"logistic","linear"),"\n",sep=""))
              cat("\n")
              
              
              tmp <- cbind(coefch[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+1:nef[k]],
                           sech[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+1:nef[k]],
                           waldch[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+1:nef[k]],
                           pwaldch[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+1:nef[k]])
              
              maxch <- apply(tmp,2,maxchar)
              if(any(c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+1:nef[k]) %in% posfix)) maxch[1] <- maxch[1]-1
              dimnames(tmp) <- list(names(coef)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+1:nef[k]],
                                    c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                      paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                      paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                      paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
              cat("\n")
              print(tmp,quote=FALSE,na.print="")
              cat("\n")
              
              # random effects
              cat("       random effects varcov matrix: \n")
              cat("\n")
              
              tmp <- coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+1:nvc[k]]
              tmp <- round(tmp,5)
              Mat.cov <- matrix(0,nrow=nea[k],ncol=nea[k])
              if(idmodel[k]==0) diag(Mat.cov) <- tmp
              if(idmodel[k]==1){
                Vect.cov <- tmp
                if(idlink[k]!=0) Vect.cov <- c(1,Vect.cov)
                Mat.cov[lower.tri(Mat.cov,diag=TRUE)] <- Vect.cov
              }
              Mat.cov[upper.tri(Mat.cov)] <- NA
              
              if(any(posfix %in% c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+1:nvc[k])))
              {
                Mat.cov <- apply(Mat.cov,2,format,digits=5,nsmall=5)
                Mat.cov[upper.tri(Mat.cov)] <- ""
                pf <- sort(intersect(c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+1:nvc[k]),posfix))
                p <- matrix(0,nea[k],nea[k])
                if(idmodel[k]==1){ 
                  if(idlink[k]==0) p[upper.tri(p,diag=TRUE)] <- c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+1:nvc[k])
                  if(idlink[k]!=0) p[upper.tri(p,diag=TRUE)] <- c(0,nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+1:nvc[k])
                }
                if(idmodel[k]==0) diag(p) <- nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+1:nvc[k]
                Mat.cov[which(t(p) %in% pf)] <- paste(Mat.cov[which(t(p) %in% pf)],"*",sep="")
              }
              
              Mat.cov.names <- NULL
              if(idmodel[k]==0){ #logistic
                Mat.cov.names <- c(Mat.cov.names,"rate.")
                Mat.cov.names <- c(Mat.cov.names,"v.")
                Mat.cov.names <- c(Mat.cov.names,"u. ")
              }
              if(idmodel[k]==1){ #linear
                Mat.cov.names <- c(Mat.cov.names,"Intercept"[which(x$idea[k]==1)])
                power.chr <- c("-2","-1","-0.5","0","0.5","1","2","3")
                for(l in 1:8){
                  if(x$ideafp[k,l]!=0){
                    Mat.cov.names <- c(Mat.cov.names,paste("FP[",power.chr[l],"]",sep=""))
                    if(x$ideafp[k,l]==2){
                      Mat.cov.names <- c(Mat.cov.names,paste("FP[",power.chr[l],".log]",sep=""))
                    }
                  }
                }
              }
              dimnames(Mat.cov) <- list(Mat.cov.names,Mat.cov.names)
                                    
              cat("\n")
              print(Mat.cov,quote=FALSE)
              cat("\n")
              
              
              
              # transformation
              if(idlink[k]!=0){
                
                if(idlink[k]==1) cat("    transformation: linear \n")
                if(idlink[k]==2) cat("    transformation: splines \n")
                cat("\n")
                
                tmp <- cbind(coefch[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+nvc[k]+1:ntr[k]],
                             sech[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+nvc[k]+1:ntr[k]],
                             waldch[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+nvc[k]+1:ntr[k]],
                             pwaldch[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+nvc[k]+1:ntr[k]])
                
                maxch <- apply(tmp,2,maxchar)
                if(any(c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+nvc[k]+1:ntr[k]) %in% posfix)) maxch[1] <- maxch[1]-1
                dimnames(tmp) <- list(names(coef)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+nvc[k]+1:ntr[k]],
                                      c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                        paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                        paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                        paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
                
                print(tmp,quote=FALSE,na.print="")
                cat("\n")
                
              }
              else{
                cat("    no transformation : identity \n")
                cat("\n")
              }
              
              
              # residual error
              cat("    error: \n")
              cat("\n")
              
              tmp <- cbind(coefch[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+nvc[k]+ntr[k]+1])
              
              maxch <- apply(tmp,2,maxchar)
              if(any(c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+nvc[k]+ntr[k]+1) %in% posfix)) maxch[1] <- maxch[1]-1
              dimnames(tmp) <- list(names(coef)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+sumPrmK+nef[k]+nvc[k]+ntr[k]+1],
                                    c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep="")))
              
              print(tmp,quote=FALSE,na.print="")
              cat("\n")
              
              
              #incrementation
              sumPrmK <- sumPrmK + nef[k]+nvc[k]+ntr[k]+1
              
            }
            
            
            if(length(posfix))
                {
                    cat(" *  coefficient fixed by the user \n \n")
                }
            
            
            return(invisible(NULL))
        }
}
