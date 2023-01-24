#' Predictions of the random-effects
#' 
#' The function computes the predicted values of the random effects given observed data
#' provided in input. 
#' 
#' @param model an object inheriting from class \code{LTSM}.
#' @param newdata data frame containing the data from which predictions are to be computed. 
#' The data frame should include at least all the covariates, 
#' the marker(s) values and the grouping structure. Names should match exactly the names 
#' of the variables in the model.
#' @param subject character specifying the name of the grouping structure.
#' If NULL (the default), the same as in the model will be used.
#' @return a matrix containing the predicted random-effects.
#' @author Tiphaine Saulnier, Cecile Proust-Lima, Viviane Philipps 
#' 
#' @export 
#' 


predRE.LTSM <- function(model, newdata, subject=NULL){
  
  if(missing(model)){ stop("The argument model should be specified")}
  if(class(model)!="LTSM") stop("The argument model must be an element inheriting from class 'LTSM'")
  if(missing(newdata)){ stop("The argument newdata should be specified")}
  if(nrow(newdata)==0) stop("newdata should not be empty")
  
  if(is.null(subject)) subject <- model$call$subject
  if(!(subject %in% colnames(newdata))) stop("Unable to find variable 'subject' in 'newdata'")
  if(!is.numeric(newdata[,subject])) stop("The argument subject must be numeric in newdata")
  nom.subject <- as.character(subject)
  
  response <- as.character(model$call$response)
  if(response[1]=="c") response <- response[-1] # if vector
  if(!all(response %in% colnames(newdata))) stop("Unable to find variable from 'response' in 'newdata'")
  range <- as.numeric(as.character(model$call$range)[-1])
  var.time <- model$call$var.time
  if(!(var.time %in% colnames(newdata))) stop("Unable to find variable 'var.time' in 'newdata'")
  fixed_s <- model$call$fixed_s
  class(fixed_s) <- "formula"
  
  data <- newdata
  
  B <- model$best
  
  # verification convergence
  if(model$conv!=1) stop("The model should have converged to predict.")
  
  
  # arguments LTSM
  # si non-renseigne : valeurs par defaut
  
  na.action <- ifelse(!is.null(model$call$na.action),model$call$na.action,1)
  
  survival <- NULL; if(!is.null(model$call$survival)) survival <- model$call$survival
  hazard <- ifelse(!is.null(model$call$hazard),model$call$hazard,"Weibull")
  hazardrange <- NULL; if(!is.null(model$call$hazardrange)) hazardrange <- model$call$hazardrange
  hazardnodes <- NULL; if(!is.null(model$call$hazardnodes)) hazardnodes <- model$call$hazardnodes
  logscale <- ifelse(!is.null(model$call$logscale),model$call$logscale,FALSE)
  startWeibull <- ifelse(!is.null(model$call$startWeibull),model$call$startWeibull,0)
  
  methInteg <- ifelse(!is.null(model$call$methInteg),model$call$methInteg,"QMC")
  nMC <- ifelse(!is.null(model$call$nMC),model$call$nMC,1000)
  maxiter <- ifelse(!is.null(model$call$maxiter),model$call$maxiter,100)
  convB <- ifelse(!is.null(model$call$convB),model$call$convB,0.0001)
  convL <- ifelse(!is.null(model$call$convL),model$call$convL,0.0001)
  convG <- ifelse(!is.null(model$call$convG),model$call$convG,0.0001)
  partialH <- NULL; if(!is.null(model$call$partialH)) partialH <- model$call$partialH
  nsim <- ifelse(!is.null(model$call$nsim),model$call$nsim,100)
  
  verbose <- ifelse(!is.null(model$call$verbose),model$call$verbose,TRUE)
  nproc <- ifelse(!is.null(model$call$nproc),model$call$nproc,1)
  clustertype <- NULL; if(!is.null(model$call$clustertype)) clustertype <- model$call$clustertype
  
  ###
  
  # print(paste("subject",subject))
  # print("response");print(paste(response))
  # print("range");print(paste(range))
  # print(paste("var.time",var.time))
  # print("fixed_s");print(paste(fixed_s))
  # print("B");print(paste(B))
  
  
  ##############################################################################
  
  #############################################
  ###   preparation des arguments Fortran   ###
  #############################################
  # strictement identique a LTSM
  
  
  ##pour acces aux attributs des formules
  
  #liste des outcomes
  nomsY <- response
  ny <- length(nomsY)
  if(!all(nomsY %in% colnames(data))) stop(paste("Data should contain the variables",paste(unique(nomsY),collapse=" ")))
  
  ## fixed_s
  var.fixed_s <- unlist(lapply(c(attr(terms(fixed_s, specials=c("factor")),"variables")),all.vars))
  if(!all(var.fixed_s %in% colnames(data))) stop(paste("Data should contain the variables",paste(unique(var.fixed_s),collapse=" ")))
  
  ##liste des variables utilisees  (sans les interactions et sans les Y)
  ttesLesVar <- c(var.time,var.fixed_s)
  
  # ## argument subset
  # form1 <- paste(c(nom.subject,nomsY,ttesLesVar),collapse="+")
  # if(!isTRUE(all.equal(as.character(cl$subset),character(0))))
  # {
  #   cc <- cl
  #   cc <- cc[c(1,which(names(cl)=="subset"))]
  #   cc[[1]] <- as.name("model.frame")
  #   cc$formula <- formula(paste("~",form1))
  #   cc$data <- data
  #   cc$na.action <- na.pass
  #   data <- eval(cc)
  # }
  # 
  # attributes(data)$terms <- NULL
  
  ## si subject est un factor
  if(is.factor(data[,nom.subject]))
  {
    data[,nom.subject] <- as.numeric(data[,nom.subject])
  }
  
  ## partie survie
  if(is.null(survival))
  {
    nbevt <- 0
    idtrunc <- 0
    nprisq <- 0
    nrisqtot <- 0
    nvarxevt <- 0
    nvarxevt2 <- 0
    typrisq <- 0
    noms.surv <- NULL
    #nom.timedepvar <- NULL
    form.commun <- ~-1
    form.cause <- ~-1
    survival <- ~-1
    nz <- 0
    zi <- 0
    minT <- 0
    maxT <- 0
  }
  else
  {
    
    stop("Survival part : not coded yet")
    
    # ## objet Surv
    # surv <- cl$survival[[2]]
    # 
    # if(length(surv)==3) #censure droite sans troncature gauche
    # {
    #     idtrunc <- 0 
    #     
    #     Tevent <- getElement(object=data,name=as.character(surv[2]))
    #     Event <- getElement(object=data,name=as.character(surv[3]))  
    #     Tentry <- rep(0,length(Tevent)) #si pas de troncature, Tentry=0
    #     
    #     noms.surv <-  c(as.character(surv[2]),as.character(surv[3]))
    #     
    #     surv <- do.call("Surv",list(time=Tevent,event=Event,type="mstate"))  
    # }
    # 
    # if(length(surv)==4) #censure droite et troncature
    # {
    #     idtrunc <- 1
    #     
    #     Tentry <- getElement(object=data,name=as.character(surv[2]))
    #     Tevent <- getElement(object=data,name=as.character(surv[3]))
    #     Event <- getElement(object=data,name=as.character(surv[4]))  
    #     
    #     noms.surv <-  c(as.character(surv[2]),as.character(surv[3]),as.character(surv[4]))   
    #     
    #     surv <- do.call("Surv",list(time=Tentry,time2=Tevent,event=Event,type="mstate"))   
    # }  
    # 
    # ## nombre d'evenement concurrents
    # nbevt <- length(attr(surv,"states"))
    # if(nbevt<1) stop("No observed event in the data")
    # 
    # 
    # ## pour la formule pour survivial, creer 3 formules : 
    # ## une pour les covariables en mixture, une pour les covariables avec effet specifique a la cause, et une pour les effets communs.  
    # form.surv <- cl$survival[3]
    # 
    # noms.form.surv <- all.vars(attr(terms(formula(paste("~",form.surv))),"variables"))
    # if(length(noms.form.surv)==0)
    # {
    #     form.cause <- ~-1
    #     form.commun <- ~-1
    #     asurv <- terms(~-1)
    # }
    # else
    # {
    #     ##creer la formula pour cause
    #     form1 <- gsub("mixture","",form.surv)       #TS: pas de classes latentes
    #     form1 <- formula(paste("~",form1))
    #     asurv1 <- terms(form1,specials="cause")  
    #     ind.cause <- attr(asurv1,"specials")$cause
    #     if(length(ind.cause))
    #     {
    #         form.cause <- paste(labels(asurv1)[ind.cause],collapse="+")
    #         form.cause <- gsub("cause","",form.cause)
    #         form.cause <- formula(paste("~",form.cause))
    #     }
    #     else
    #     {
    #         form.cause <- ~-1 
    #     }
    #     
    #     
    #     ## creer la formule pour ni cause ni mixture
    #     asurv <- terms(formula(paste("~",form.surv)),specials=c("cause"))
    #     ind.commun <- setdiff(1:length(labels(asurv)),unlist(attr(asurv,"specials")))
    #     if(length(ind.commun))
    #     {
    #         form.commun <- paste(labels(asurv)[ind.commun],collapse="+")
    #         form.commun <- gsub("cause","",form.commun)   # si X1:cause(X2)
    #         form.commun <- formula(paste("~",form.commun))  
    #     }
    #     else
    #     {
    #         form.commun <- ~-1 
    #     }
    # }
    #
    # ##verifier si toutes les variables sont dans data
    # varSurv <- unique(all.vars(terms(survival)))
    # if(!all(varSurv %in% colnames(data))) stop(paste("Data should contain the variables",paste(varSurv,collapse=" ")))
    # ttesLesVar <- unique(c(ttesLesVar,varSurv))
  }
  
  
  ##subset de data avec les variables utilisees
  newdata <- data[,c(nom.subject,nomsY,noms.surv,ttesLesVar),drop=FALSE]
  
  dataSurv <- NULL
  if((nbevt>0))
  {
    dataSurv <- data.frame(getElement(object=data,name=nom.subject),Tentry,Tevent,Event)
  }
  
  
  ##un data frame par outcome et creation Y0
  dataY <- paste("data",nomsY,sep=".")
  Y0 <- NULL
  IND <- NULL
  outcome <- NULL
  data0 <- NULL
  nayk <- vector("list",ny)
  for (k in 1:ny)
  {
    dtemp <- newdata[,c(nom.subject,nomsY[k],noms.surv,ttesLesVar)]
    ##enlever les NA
    linesNA <- apply(dtemp,2,function(v) which(is.na(v)))
    linesNA <- unique(unlist(linesNA))
    if(length(linesNA)) nayk[[k]] <- linesNA
    if(na.action==1 & length(linesNA)>0) dtemp <- dtemp[-linesNA,]
    if(na.action==2 & length(linesNA)>0) stop("Data contains missing values")
    assign(dataY[k],dtemp)
    Y0 <- c(Y0, dtemp[,nomsY[k]])
    IND <- c(IND, dtemp[,nom.subject])
    outcome <- c(outcome,rep(nomsY[k],nrow(dtemp)))
    data0 <- rbind(data0, dtemp[,setdiff(colnames(dtemp),nomsY[k]),drop=FALSE])   #dataset sans NA avec les covariables utilisees; obs ordonnees par outcome
  }
  
  
  ##creation de X0 (ttes les var + interactions)
  
  # vecteur des temps
  Xtime <- model.matrix(formula(paste("~",var.time)), data=data0)[,2]
  
  # matrice varexp pour s
  if(length(var.fixed_s)>0){ # verif varexp time independent
    for(v in 1:length(var.fixed_s)){ #verif: time independent
      tmp <- unique(na.omit(data0[,c(subject,var.fixed_s[v])]))  #dataframe 2 colonnes : subject et var v, en supprimant les lignes doublons
      if(nrow(tmp) != length(unique(IND))) #var v dependante du temps
        stop(paste("In fixed_s, variable ",var.fixed_s[v]," seems to be time dependant, it can't be"))
    }  
  }
  Xfixed_s <- model.matrix(fixed_s, data=data0)[,-1] # sans intercept
  
  # matrices de survie
  Xsurv <- model.matrix(form.commun,data=data0)  #+
  Xsurvcause <- model.matrix(form.cause,data=data0)  #+
  
  z.fixed_s <- strsplit(colnames(Xfixed_s),split=":",fixed=TRUE)
  z.fixed_s <- lapply(z.fixed_s,sort)
  
  if(form.commun != ~-1)
  {
    z.surv <- strsplit(colnames(Xsurv),split=":",fixed=TRUE)
    z.surv <- lapply(z.surv,sort)
  }
  else
  {
    z.surv <- list() 
  }
  
  if(form.cause != ~-1)
  {
    z.survcause <- strsplit(colnames(Xsurvcause),split=":",fixed=TRUE)
    z.survcause <- lapply(z.survcause,sort)
  }
  else
  {
    z.survcause <- list() 
  }
  
  X0 <- cbind(Xtime, Xfixed_s, Xsurv, Xsurvcause)
  
  colnames(X0)[1] <- var.time
  nom.unique <- unique(colnames(X0))
  X0 <- X0[,nom.unique,drop=FALSE]
  X0 <- as.matrix(X0)
  ##X0 fini
  
  
  ## transfo outcomes
  # argument range
  if(length(range) != 2*ny) stop("Length of vector range is not correct.")
  for(k in 1:ny){
    rg_thq <- range[2*(k-1)+1:2]
    rg_obs <- range(get(dataY[k])[,nomsY[k]])
    if(rg_obs[1]<rg_thq[1] | rg_obs[2]>rg_thq[2]) stop("The range specified do not cover the entire range of the data")
  }
  #remplir zitr
  zitr <- matrix(range,2,ny)
  
  
  ##ordonner les mesures par individu
  matYX <- cbind(IND,Y0,outcome,X0)
  matYXord <- matYX[order(IND),]
  Y0 <- as.numeric(matYXord[,2])
  X0 <- apply(matYXord[,-c(1:3),drop=FALSE],2,as.numeric)
  IND <- matYXord[,1]
  outcome <- matYXord[,3]
  
  if(survival != ~-1){
    
    stop("Survival part : not coded yet")
    
    # dataSurv <- dataSurv[which(dataSurv[,1] %in% IND),]
    # dataSurv <- dataSurv[order(dataSurv[,1]),]
    # nmes <- as.vector(table(dataSurv[,1]))
    # data.surv <- apply(dataSurv[cumsum(nmes),-1],2,as.numeric)
    # tsurv0 <- data.surv[,1]
    # tsurv <- data.surv[,2]
    # devt <- data.surv[,3]
    # 
    # 
    # ## test de hazard
    # arghaz <- hazard
    # hazard <- rep(hazard,length.out=nbevt)
    # if(any(hazard %in% c("splines","Splines")))
    # {
    #   hazard[which(hazard %in% c("splines","Splines"))] <- "5-quant-splines"
    # }
    # if(any(hazard %in% c("piecewise","Piecewise")))
    # {
    #   hazard[which(hazard %in% c("piecewise","Piecewise"))] <- "5-quant-piecewise"
    # }
    # 
    # haz13 <- strsplit(hazard[which(!(hazard=="Weibull"))],"-")
    # if(any(sapply(haz13,length)!=3)) stop("Invalid argument hazard")
    # 
    # nz <- rep(2,nbevt)
    # locnodes <- NULL
    # typrisq <- rep(2,nbevt)
    # nprisq <- rep(2,nbevt)
    # 
    # nznodes <- 0 #longueur de hazardnodes
    # ii <- 0
    # if(any(hazard!="Weibull"))
    # {
    #   
    #   for (i in 1:nbevt)
    #   {
    #     if(hazard[i]=="Weibull") next;
    #     
    #     ii <- ii+1
    #     
    #     nz[i] <- as.numeric(haz13[[ii]][1])
    #     if(nz[i]<3) stop("At least 3 nodes are required")
    #     typrisq[i] <- ifelse(haz13[[ii]][3] %in% c("splines","Splines"),3,1)
    #     nprisq[i] <- ifelse(haz13[[ii]][3] %in% c("splines","Splines"),nz[i]+2,nz[i]-1)
    #     locnodes <- c(locnodes, haz13[[ii]][2])
    #     if(!(haz13[[ii]][3] %in% c("splines","Splines","piecewise","Piecewise"))) stop("Invalid argument hazard")
    #     
    #     if((haz13[[ii]][2]=="manual"))
    #     {
    #       nznodes <- nznodes + nz[i]-2
    #     }
    #     
    #     if(!all(locnodes %in% c("equi","quant","manual"))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")
    #   }
    #   
    #   if(!is.null(hazardnodes))
    #   {
    #     if(!any(locnodes == "manual"))  stop("hazardnodes should be NULL if the nodes are not chosen manually")
    #     
    #     if(length(hazardnodes) != nznodes) stop(paste("Vector hazardnodes should be of length",nznodes))
    #   }
    # }
    # else
    # {
    #   if(!is.null(hazardnodes)) stop("hazardnodes should be NULL if Weibull baseline risk functions are chosen")
    # }
    # 
    # 
    # if(nbevt>1 & length(arghaz)==1 & nznodes>0)
    # {
    #   hazardnodes <- rep(hazardnodes,length.out=nznodes*nbevt)
    # }
    # 
    # nrisqtot <- sum(nprisq)
    # 
    # zi <- matrix(0,nrow=max(nz),ncol=nbevt)
    # nb <- 0
    # 
    # minT1 <- 0
    # maxT1 <- max(tsurv)
    # tsurvevt <- tsurv
    # 
    # if(idtrunc==1)
    # {
    #   minT1 <- min(tsurv,tsurv0)
    #   maxT1 <- max(tsurv,tsurv0)
    # }
    # 
    # ## arrondir
    # minT2 <- round(minT1,3)
    # if(minT1<minT2) minT2 <- minT2-0.001
    # minT <- minT2
    # 
    # maxT2 <- round(maxT1,3)
    # if(maxT1>maxT2) maxT2 <- maxT2+0.001
    # maxT <- maxT2
    # 
    # if(length(hazardrange)){
    #   if(hazardrange[1]>minT) stop(paste("hazardrange[1] should be <=",minT))
    #   if(hazardrange[2]>maxT) stop(paste("hazardrange[2] should be >=",maxT))
    #   minT <- hazardrange[1]
    #   maxT <- hazardrange[2]
    # }
    # 
    # startWeib <- rep(0,nbevt)
    # startWeib[which(typrisq==2)] <- rep(startWeibull, length.out=length(which(typrisq==2)))
    # ii <- 0
    # for(i in 1:nbevt)
    # {
    #   if(typrisq[i]==2)
    #   {
    #     if(minT < startWeib[i]) stop("Some entry or event times are bellow startWeibull")
    #     zi[1:2,i] <- c(startWeib[i],maxT)
    #   }
    #   else
    #   {
    #     ii <- ii+1
    #     
    #     if(locnodes[ii]=="manual")
    #     {
    #       zi[1:nz[i],i] <- c(minT,hazardnodes[nb+1:(nz[i]-2)],maxT)
    #       nb <- nb + nz[i]-2
    #     }
    #     if(locnodes[ii]=="equi")
    #     {
    #       zi[1:nz[i],i] <- seq(minT,maxT,length.out=nz[i])
    #     }
    #     if(locnodes[ii]=="quant")
    #     {
    #       pi <- c(1:(nz[i]-2))/(nz[i]-1)
    #       qi <- quantile(tsurvevt,prob=pi)
    #       zi[1,i] <- minT
    #       zi[2:(nz[i]-1),i] <- qi
    #       zi[nz[i],i] <- maxT
    #     }
    #   }
    # }
  }
  else{
    tsurv0 <- rep(0,length(unique(IND)))
    tsurv <- rep(0,length(unique(IND))) 
    devt <- rep(0,length(unique(IND)))
  }
  
  
  ##parametres pour Fortran
  ns <- length(unique(IND))
  nv <- dim(X0)[2]
  nobs <- length(Y0)
  logspecif <- as.numeric(logscale)
  #loglik <- 0
  ni <- 0
  istop <- 0
  gconv <- rep(0,3)
  resid_m <- rep(0,nobs)
  resid_ss <- rep(0,nobs)
  Yobs <- rep(0,nobs)
  time <- seq(minT,maxT,length.out=nsim)
  risq_est <- matrix(0,nrow=nsim,ncol=nbevt)
  risqcum_est <- matrix(0,nrow=nsim,ncol=nbevt)
  #predRE_Y <- rep(0,ns*1*ny)  #+
  rlindiv <- rep(0,ns)
  marker <- rep(0,nsim*ny)
  transfY <- rep(0,nsim*ny)
  
  
  ##nmes
  nmes <- matrix(0,ns,ny)
  for (k in 1:ny)
  {
    INDpresents <- which(unique(IND) %in% get(dataY[k])[,nom.subject])
    nmes[INDpresents,k] <- as.vector(table(get(dataY[k])[,nom.subject]))
  }
  maxmes <- max(apply(nmes,1,sum))
  
  
  ##remplir idprob, etc
  z.X0 <- strsplit(nom.unique,split=":",fixed=TRUE)
  z.X0 <- lapply(z.X0,sort)
  
  idg <- (z.X0 %in% z.fixed_s) + 0
  
  idsurv <- z.X0 %in% z.surv + 2*(z.X0 %in% z.survcause)
  idsurv[1] <- 0 # 0 pour l'intercept
  
  ## nb coef ppour survie
  nvarxevt <- sum(idsurv==1) + nbevt*sum(idsurv==2)
  
  
  ## nb parametres d'association
  nasso <- (1+2*ny)*nbevt
  
  
  ## prm partie long
  nef.s <- sum(idg==1) #time shift
  nvc.s <- 1
  nef.B <- ny #rate
  nmu.v <- ny #location
  nvc.v <- ny
  nvc.u <- ny #intercept aleatoire
  nvc.err <- ny #erreur
  
  ##nombre total de parametres
  NPM <- nrisqtot + nvarxevt + nasso + 
    nef.s + nvc.s + nef.B + nmu.v + nvc.v + nvc.u + nvc.err
  
  ##nombre effets aleatoires
  nea <- nvc.s + nvc.v + nvc.u
  
  
  # V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres
  
  
  ###valeurs initiales
  # if(!(missing(B)))
  # {
  #   if(!is.vector(B)) stop("B should be a vector")
  #   
    if (length(B)==NPM) b <- B
    else stop(paste("Vector B should be of length",NPM))
  #   
  # }
  # else ## B missing
  # {
  #   b <- rep(0,NPM)
  #   
  #   if(nbevt>0){  #TS : si Weibull et prms au carre pr positivite -> valeurs par defaut = 1
  #     if(any(hazard!="Weibull")==FALSE & isFALSE(logscale)==TRUE){
  #       for(i in 1:nbevt){
  #         if(typrisq[i]==2){
  #           b[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- 1
  #           if(idtrunc==1 & any(Tentry==0)) #TS : si entree retardee et au moins un temps vaut 0
  #             b[sum(nprisq[1:i])-nprisq[i]+nprisq[i]] <- 1.25 #sinon pblm lors de recherche prms, risq instant tend vers infini qd w2 < 1 et t=0
  #         }
  #       }
  #     }
  #     if(any(hazard %in% c("splines","Splines")) & idtrunc==1 & any(Tentry==0)){
  #       for(i in 1:nbevt){
  #         if(typrisq[i]==3){
  #           b[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- 10^-7 #sinon pgrm crash
  #         }
  #       }
  #     }
  #   }
  #   
  #   #1 for variances
  #   b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s] <- 1
  #   b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + 1:nvc.v] <- rep(1,nmu.v)
  #   b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + 1:nvc.u] <- rep(1,nvc.u)
  #   b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + nvc.u + 1:nvc.err] <- rep(1,nvc.err)
  #   
  # }
  
  
  # ##------------------------------------------
  # ##------nom au vecteur best
  # ##--------------------------------------------
  # 
  # nom.X0 <- colnames(X0)
  # nom.X0[nom.X0=="(Intercept)"] <- "intercept"
  # 
  # if(nbevt>0)
  # {
  #   ##prm fct de risque
  #   if(isTRUE(logscale))
  #   {
  #     for(i in 1:nbevt)
  #     {
  #       nom1 <- rep(paste("event",i,sep=""),nprisq[i])
  #       if(typrisq[i]==2)
  #       {
  #         names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:2]," log(Weibull",1:2,")",sep="")
  #       }
  #       if(typrisq[i]==1)
  #       {
  #         names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(piecewise",1:(nz[i]-1),")",sep="")
  #       }
  #       if(typrisq[i]==3)
  #       {
  #         names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(splines",1:(nz[i]+2),")",sep="")
  #       }
  #     }
  #   }
  #   else
  #   {
  #     for(i in 1:nbevt)
  #     {
  #       nom1 <- rep(paste("event",i,sep=""),nprisq[i])
  #       if(typrisq[i]==2)
  #       {
  #         names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:2]," +/-sqrt(Weibull",1:2,")",sep="")
  #       }
  #       if(typrisq[i]==1)
  #       {
  #         names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(piecewise",1:(nz[i]-1),")",sep="")
  #       }
  #       if(typrisq[i]==3)
  #       {
  #         names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(splines",1:(nz[i]+2),")",sep="")
  #       }
  #     }
  #   }
  #   
  #   
  #   ##prm covariables survival
  #   nom1 <- NULL
  #   for(j in 1:nv)
  #   {
  #     if(idsurv[j]==1) #X
  #     {
  #       if(idtdv[j]==1)
  #       {
  #         nom1 <- c(nom1,paste("I(t>",nom.timedepvar,")",sep=""))
  #       }
  #       else
  #       {
  #         nom1 <- c(nom1,nom.X0[j])
  #       }
  #     }
  #     
  #     if(idsurv[j]==2) #cause(X)
  #     {
  #       if(idtdv[j]==1)
  #       {
  #         nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") event",1:nbevt,sep=""))
  #       }
  #       else
  #       {
  #         nom1 <- c(nom1,paste(nom.X0[j],paste("event",1:nbevt,sep="")))
  #       }
  #     }
  #     
  #   }
  #   
  #   if(nvarxevt>0) names(b)[nrisqtot+1:nvarxevt] <- nom1
  #   
  #   
  #   for(i in 1:nbevt)
  #     names(b)[nrisqtot+nvarxevt+(nbevt-1)*nea+1:nea] <- paste("event",i," asso",1:nea,sep="")
  #   
  # }
  # 
  # 
  # # partie long
  # names(b)[nrisqtot+nvarxevt+nasso+1:nef.s] <- nom.X0[idg!=0]
  # names(b)[nrisqtot+nvarxevt+nasso+nef.s+1] <- "var.s"
  # names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+1:nef.B] <- paste("rate",1:ny,sep="")
  # names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+1:nmu.v] <- paste("mu.v",1:ny,sep="")
  # names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+1:nvc.v] <- paste("var.v",1:ny,sep="")
  # names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+1:nvc.u] <- paste("var.u",1:ny,sep="")
  # names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+nvc.u+1:nvc.err] <- paste("std.err",1:ny,sep="")
  
  namesb <- names(b)
  
  
  # ## prm fixes
  # fix <- rep(0,NPM)
  # if(length(posfix))
  # {
  #   if(any(!(posfix %in% 1:NPM))) stop("Indexes in posfix are not correct")
  #   
  #   fix[posfix] <- 1
  # }
  # if(length(posfix)==NPM) stop("No parameter to estimate")
  # 
  # if(!all(partialH %in% 1:NPM)) stop(paste("partialH should contain indices between 1 and",NPM))
  
  # RE : var -> std
  b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s] <- sqrt(b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s])
  b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + 1:nvc.v] <- sqrt(b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + 1:nvc.v])
  b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + 1:nvc.u] <- sqrt(b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + 1:nvc.u])
  #b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + nvc.u + 1:nvc.err] <- sqrt(b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + nvc.u + 1:nvc.err])
  
  # # fixed prms
  # nfix <- sum(fix)
  # bfix <- 0
  # if(nfix>0)
  # {
  #   bfix <- b[which(fix==1)]
  #   b <- b[which(fix==0)]
  #   NPM <- NPM-nfix
  # }
  
  
  ##############################################################################
  
  ##predictions
  predRE <- matrix(NA,nrow=ns,ncol=nea)  # 1 pr s, 1 / mq pr u, 1 / mq pr v
  predRE <- data.frame(unique(IND),rep(NA,ns),rep(NA,ns),predRE)
  colnames(predRE) <- c(nom.subject,"log.dens.value","conv","re.s",paste("re.v",1:ny),paste("re.u",1:ny))
  
  ##############################################################################
  
  ###############
  ###   MLA   ###
  ###############
  
  
  # boucle sur les sujets
  # for i
  #      mla
  #      sauvegarde pred
  
  # valeurs initiales
  re <- rep(1, nea)
  nmescur <- 0
  
  for(i in 1:ns){
    
    #print("#########################")
    #print(paste("SUBJECT ",i))
    
    res <- mla(b = re, m = nea, fn = predREfct,
               clustertype=clustertype,.packages=NULL,
               epsa=convB,epsb=convL,epsd=convG,
               digits=8,print.info=verbose,blinding=FALSE,
               multipleTry=25,file="",partialH=partialH,
               nproc=nproc,maxiter=maxiter,minimize=FALSE,
               best0=b,
               Y0=Y0,X0=X0,idg0=idg,
               ny0=ny,ns0=ns,nv0=nv,
               nobs0=nobs,nmes0=nmes,npm0=NPM,nea0=nea,
               zitr0=zitr,
               ind0=i,nmescur0=nmescur)
    
    # argument nmescur, ...
    
    
    predRE[i,2] <- res$fn.value # fct density value
    predRE[i,3] <- res$istop # conv
    predRE[i,-c(1:3)] <- res$b # predRE
    
    
    nmescur <- nmescur + sum(nmes[i,])
    
    # print(paste("dens value ",res$fn.value))
    # print(paste("conv ",res$istop))
    # print("re ");print(paste(res$b))
    # print(paste("nmescur_i ",nmescur))
    
  }
  
  
  ##############################################################################
  
  res <-list(ns=ns, predRE=predRE, best=B)
  
  return(res)
  
} 

predREfct <- function(re0,best0,Y0,X0,idg0,
                   ny0,ns0,nv0,nobs0,nmes0,npm0,nea0,zitr0,
                   ind0,nmescur0)
{
  res <- 0
  .Fortran(C_dens,as.double(re0),as.double(best0),as.double(Y0),as.double(X0),as.integer(idg0),as.integer(ny0),as.integer(ns0),as.integer(nv0),as.integer(nobs0),as.integer(nmes0),as.integer(npm0),as.integer(nea0),as.double(zitr0),as.integer(ind0),as.integer(nmescur0),dens_res=as.double(res))$dens_res
}
