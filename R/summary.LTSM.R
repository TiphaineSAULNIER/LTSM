#' @export
#'
summary.LTSM <- function(object,...)
{
    x <- object

    cat("Latent time shift model fitted by maximum likelihood method \n")
    
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

    cat(paste("     Number of observations:", x$N[13]),"\n")
    cat(paste("     Number of parameters:", length(x$best))," \n")
    if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")


    #nbevt <- x$N[14]
    # if(nbevt>0)
    # {
    #     nprisq <- rep(NA,nbevt)
    #     nrisq <- rep(NA,nbevt)
    #     for(ke in 1:nbevt)
    #     {
    #         if(x$typrisq[ke]==1) nprisq[ke] <- x$nz[ke]-1
    #         if(x$typrisq[ke]==2) nprisq[ke] <- 2
    #         if(x$typrisq[ke]==3) nprisq[ke] <- x$nz[ke]+2
    # 
    #         
    #         cat(paste("     Event",ke,": \n"))
    #         cat(paste("        Number of events: ", x$nevent[ke],"\n",sep=""))
    #         
    #         if (x$typrisq[ke]==2)
    #         {
    #             cat("        Weibull baseline risk function \n")
    #         }
    #         if (x$typrisq[ke]==1)
    #         {
    #             cat("        Piecewise constant baseline risk function with nodes \n")
    #             cat("        ",x$hazardnodes[1:x$nz[ke],ke]," \n")
    #         }
    #         if (x$typrisq[ke]==3)
    #         {
    #             cat("        M-splines constant baseline risk function with nodes \n")
    #             cat("        ",x$hazardnodes[1:x$nz[ke],ke]," \n")
    #         }
    #                     
    #     }
    # }

    
    # ny <- x$N[12]
    # ntr <- rep(NA,ny)
    # numSPL <- 0
    # cat("     Link functions: ")
    # for (yk in 1:ny)
    #     {
    #         if (x$linktype[yk]==0)
    #             {
    #                 ntr[yk] <- 2
    #                 if (yk>1) cat("                     ")
    #                 cat("Linear for",x$Names$Ynames[yk]," \n")
    #             }
    #         if (x$linktype[yk]==1)
    #             {
    #                 ntr[yk] <- 4
    #                 if (yk>1) cat("                     ")
    #                 cat("Standardised Beta CdF for",x$Names$Ynames[yk]," \n")
    #             }
    #         if (x$linktype[yk]==2) 
    #             {
    #                 numSPL <- numSPL+1
    #                 ntr[yk] <- x$nbnodes[numSPL]+2
    #                 if (yk>1) cat("                     ")
    #                 cat("Quadratic I-splines with nodes", x$linknodes[1:x$nbnodes[numSPL],yk],"for",x$Names$Ynames[yk], "\n")
    #             }
    #         if (x$linktype[yk]==3) 
    #             {
    #                 ntr[yk] <- x$nbmod[yk]-1
    #                 if (yk>1) cat("                     ")
    #                 cat("Thresholds for",x$Names$Ynames[yk], "\n")
    #             }
    #     }

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

            nprob <- x$N[1]
            nrisqtot <- x$N[2]
            nvarxevt <- x$N[3]
            nasso <- x$N[4]
            nef.s <- x$N[5]
            nvc.s <- x$N[6]
            nef.B <- x$N[7]
            nmu.v <- x$N[8]
            nvc.v <- x$N[9]
            nvc.u <- x$N[10]
            nvc.err <- x$N[11]
            ny <- x$N[12]
            nobs <- x$N[13]
            nbevt <- x$N[14]
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

            coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+nvc.u+1:nvc.err] <- abs(coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+nvc.u+1:nvc.err]) # std.err
            
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

            # if(nbevt>0)
            # {
            #     cat("Survival model:\n" )
            #     cat("\n")
            #     cat("Parameters in the proportional hazard model:\n" )
            #     cat("\n")
            #     
            #     tmp <- cbind(coefch[c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso)],
            #                  sech[c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso)],
            #                  waldch[c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso)],
            #                  pwaldch[c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso)])
            #     maxch <- apply(tmp,2,maxchar)
            #     if(any(c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso) %in% posfix)) maxch[1] <- maxch[1]-1
            #     dimnames(tmp) <- list(names(coef)[1:(nrisqtot+nvarxevt+nasso)],
            #                           c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
            #                             paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
            #                             paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
            #                             paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
            #     cat("\n")
            #     print(tmp,quote=FALSE,na.print="")
            #     cat("\n")
            # }

            ##############################

            cat("Time shift:\n" )
            cat(" \n")
            
            tmp <- NULL
            if (nef.s>0)
            {
              tmp <- cbind(round(coef[nrisqtot+nvarxevt+nasso+1:nef.s],5),round(se[nrisqtot+nvarxevt+nasso+1:nef.s],5),round(wald[nrisqtot+nvarxevt+nasso+1:nef.s],3),round(pwald[nrisqtot+nvarxevt+nasso+1:nef.s],5))
              dimnames(tmp) <- list(names(coef)[nrisqtot+nvarxevt+nasso+1:nef.s], c("coef", "Se", "Wald", "p-value"))
            }
            
            rownames.tmp <- rownames(tmp)
            
            if((nef.s>0) & any(c(nrisqtot+nvarxevt+nasso+1:nef.s) %in% posfix)) # if 1 prm fixed
            {      
              col1 <- rep(NA,length(tmp[,1]))
              col1[which(!is.na(tmp[,1]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmp[,1]))),nsmall=5,scientific=FALSE)
              col2 <- rep(NA,length(tmp[,2]))
              col2[which(!is.na(tmp[,2]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmp[,2]))),nsmall=5,scientific=FALSE)
              col3 <- rep(NA,length(tmp[,3]))
              col3[which(!is.na(tmp[,3]))] <- format(as.numeric(sprintf("%.3f",na.omit(tmp[,3]))),nsmall=3,scientific=FALSE)
              col4 <- rep(NA,length(tmp[,4]))
              col4[which(!is.na(tmp[,4]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmp[,4]))),nsmall=5,scientific=FALSE)
              
              pf <- sort(intersect(c(nrisqtot+nvarxevt+1:(nef.s)),posfix))
              p <- rep(0,length(tmp[,1]))
              p[which(rownames(tmp) %in% c(x$Names$Xnames,x$Names$Ynames[-ny]))] <- c(nrisqtot+nvarxevt+1:(nef.s))
              col1[which(p %in% pf)] <- paste(col1[which(p %in% pf)],"*",sep="")
              col2[which(p %in% pf)] <- NA
              col3[which(p %in% pf)] <- NA
              col4[which(p %in% pf)] <- NA
              
              tmp <- cbind(col1,col2,col3,col4)
              rownames(tmp) <- rownames.tmp
              maxch <- apply(tmp,2,maxchar)
              maxch[1] <- maxch[1]-1
              
              colnames(tmp) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                 paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                 paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                 paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep=""))
            }
            
            colnames.tmp <- colnames(tmp)
            
            tmp <- rbind(tmp,
                         c(round(coef[nrisqtot+nvarxevt+nasso+nef.s+1:nvc.s],5),NA,NA,NA))
            
            dimnames(tmp) <- list(c(rownames.tmp,"Residual variance of the time shift"), colnames.tmp)
            prmatrix(round(tmp,5),na.print="")
            
            cat("\n")
            
            ##########
            
            cat("Structural Sigmoïd Model:\n" )
            cat("\n")
            
            # rate of progression
            
            tmpR <- cbind(round(coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+1:nef.B],5),round(se[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+1:nef.B],5),round(wald[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+1:nef.B],3),round(pwald[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+1:nef.B],5))
            dimnames(tmpR) <- list(paste("Mean rate ",x$Names$Ynames,sep=""), c("coef", "Se", "Wald", "p-value"))
            
            rownames.tmpR <- rownames(tmpR)
            
            if(any(c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+1:nef.B) %in% posfix)) # if 1 prm fixed
            {      
              col1 <- rep(NA,length(tmpR[,1]))
              col1[which(!is.na(tmpR[,1]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmpR[,1]))),nsmall=5,scientific=FALSE)
              col2 <- rep(NA,length(tmpR[,2]))
              col2[which(!is.na(tmpR[,2]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmpR[,2]))),nsmall=5,scientific=FALSE)
              col3 <- rep(NA,length(tmpR[,3]))
              col3[which(!is.na(tmpR[,3]))] <- format(as.numeric(sprintf("%.3f",na.omit(tmpR[,3]))),nsmall=3,scientific=FALSE)
              col4 <- rep(NA,length(tmpR[,4]))
              col4[which(!is.na(tmpR[,4]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmpR[,4]))),nsmall=5,scientific=FALSE)
              
              pf <- sort(intersect(c(nrisqtot+nvarxevt+nef.s+nvc.s+1:nef.B),posfix))
              p <- rep(0,length(tmpR[,1]))
              p[which(rownames(tmpR) %in% c(x$Names$Xnames,x$Names$Ynames[-ny]))] <- c(nrisqtot+nvarxevt+nef.s+nvc.s+1:nef.B)
              col1[which(p %in% pf)] <- paste(col1[which(p %in% pf)],"*",sep="")
              col2[which(p %in% pf)] <- NA
              col3[which(p %in% pf)] <- NA
              col4[which(p %in% pf)] <- NA
              
              tmpR <- cbind(col1,col2,col3,col4)
              rownames(tmpR) <- rownames.tmpR
              maxch <- apply(tmpR,2,maxchar)
              maxch[1] <- maxch[1]-1
              
              colnames(tmpR) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                 paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                 paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                 paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep=""))
            }
            
            
            # location shift
            
            tmpL <- cbind(round(coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+1:nmu.v],5),round(se[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+1:nmu.v],5),round(wald[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+1:nmu.v],3),round(pwald[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+1:nmu.v],5))
            dimnames(tmpL) <- list(paste("Mean location shift ",x$Names$Ynames,sep=""), c("coef", "Se", "Wald", "p-value"))
            
            rownames.tmpL <- rownames(tmpL)
            
            if(any(c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+1:nmu.v) %in% posfix)) # if 1 prm fixed
            {      
              col1 <- rep(NA,length(tmpL[,1]))
              col1[which(!is.na(tmpL[,1]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmpL[,1]))),nsmall=5,scientific=FALSE)
              col2 <- rep(NA,length(tmpL[,2]))
              col2[which(!is.na(tmpL[,2]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmpL[,2]))),nsmall=5,scientific=FALSE)
              col3 <- rep(NA,length(tmpL[,3]))
              col3[which(!is.na(tmpL[,3]))] <- format(as.numeric(sprintf("%.3f",na.omit(tmpL[,3]))),nsmall=3,scientific=FALSE)
              col4 <- rep(NA,length(tmpL[,4]))
              col4[which(!is.na(tmpL[,4]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmpL[,4]))),nsmall=5,scientific=FALSE)
              
              pf <- sort(intersect(c(nrisqtot+nvarxevt+nef.s+nvc.s+nef.B+1:nmu.v),posfix))
              p <- rep(0,length(tmpL[,1]))
              p[which(rownames(tmpL) %in% c(x$Names$Xnames,x$Names$Ynames[-ny]))] <- c(nrisqtot+nvarxevt+nef.s+nvc.s+nef.B+1:nmu.v)
              col1[which(p %in% pf)] <- paste(col1[which(p %in% pf)],"*",sep="")
              col2[which(p %in% pf)] <- NA
              col3[which(p %in% pf)] <- NA
              col4[which(p %in% pf)] <- NA
              
              tmpL <- cbind(col1,col2,col3,col4)
              rownames(tmpL) <- rownames.tmpL
              maxch <- apply(tmpL,2,maxchar)
              maxch[1] <- maxch[1]-1
              
              colnames(tmpL) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                  paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                  paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                  paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep=""))
            }
            
            colnames.tmpL <- colnames(tmpL)
            
            for(k in 1:nvc.v)
              tmpL <- rbind(tmpL,
                            c(round(coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+k],5),NA,NA,NA))
            
            dimnames(tmpL) <- list(c(rownames.tmpL,paste("Residual variance of the location shift ",x$Names$Ynames,sep="")), colnames.tmpL)
            
            
            # combine tmpR and tmpL
            
            tmp <- rbind(tmpR,tmpL)
            
            maxch <- apply(tmp,2,maxchar)
            if(any(c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+1:(nef.B+nmu.v+nvc.v)) %in% posfix))
            {
              maxch[grep("*",tmp[1,])] <- maxch[grep("*",tmp[1,])]-1
            }
            
            prmatrix(round(tmp,5),na.print="")
            cat("\n")
            
            ##########
            
            cat("Measurement Model:\n" )
            cat("\n")
            
            # random intercept
            tmpI <- cbind(coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+1:nvc.u])
            rownames(tmpI) <- paste("Residual variance of the random intercept ",x$Names$Ynames,sep="")
            colnames(tmpI) <- "coef"
            
            # residual error
            tmpE <- cbind(coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+nvc.u+1:nvc.err])
            rownames(tmpE) <- paste("Standard deviation of the residual error ",x$Names$Ynames,sep="")
            colnames(tmpE) <- "coef"
            
            
            # combine tmpI and tmpE
            
            tmp <- rbind(tmpI,tmpE)
            
            maxch <- apply(tmp,2,maxchar)
            if(any(c(nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+1:(nvc.u+nvc.err)) %in% posfix))
            {
              maxch[grep("*",tmp[1,])] <- maxch[grep("*",tmp[1,])]-1
            }
            
            prmatrix(round(tmp,5),na.print="")
            cat("\n")
            
            ##############################
            
            # cat("\n")
            # cat("Random effects (recap):\n" )
            # cat("\n")
            # 
            # Mat.cov <- diag(c(coef[nrisqtot+nvarxevt+nasso+nef.s+1:nvc.s],
            #                   coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+1:nvc.v],
            #                   coef[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+1:nvc.u]))
            # Mat.cov[lower.tri(Mat.cov)] <- 0
            # Mat.cov[upper.tri(Mat.cov)] <- NA
            # 
            # colnames(Mat.cov) <- c("s",paste("v.",x$Names$Ynames,sep=""),paste("u.",x$Names$Ynames,sep=""))
            # rownames(Mat.cov) <- colnames(Mat.cov)
            # 
            # prmatrix(round(Mat.cov,5),na.print="")
            # 
            # cat("\n")
            
            ######################
            
            if(length(posfix))
                {
                    cat(" *  coefficient fixed by the user \n \n")
                }
            
            
            return(invisible(NULL))
        }
}
