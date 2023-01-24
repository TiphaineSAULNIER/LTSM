#' Estimation of latent time shift model for multivariate longitudinal 
#' outcomes (and time-to-event data) with a time-shifting process.
#' 
#' This function fits XXX.
#' 
#' 
#' 
#' A. GENERALIZED LOGISTIC MODELS
#' 
#' XXX
#' multivariate version : one sub-model per outcome
#' 
#' 
#' 
#' B. THE SURVIVAL MODEL
#' 
#' association btw longitudinal and survival data : shared random effect (?)
#' 
#' 
#' 
#' C. THE VECTOR OF PARAMETERS B
#' 
#' The parameters in the vector of initial values \code{B} or in the vector of
#' maximum likelihood estimates \code{best} are included in the following
#' order:
#' (1) parameters for the baseline risk function: 2 parameters for each Weibull, 
#' nz-1 for each piecewise constant risk and nz+2 for each splines risk. In the 
#' presence of competing events, the number of parameters should be adapted to 
#' the number of causes of event;
#' (2) for all covariates in survival, one parameter 
#' is required. Covariates parameters should be included in the same order as in survival.
#' In the presence of cause-specific effects, the number of parameters should be
#' multiplied by the number of causes;
#' (3) parameter(s) of association between the longitudinal 
#' and the survival process: one parameter per random effect
#' and per cause of event is required;
#' (4) parameters for the time-shift s_i: one parameter per covariate in \code{fixed_s} as a fixed effect
#' (no intercept), and one variance parameter as the random intercept;
#' (5) the outcome-specific rate B_k: one parameter per response-outcome;
#' (6) the outcome-specific random location v_ik: one mean parameter per outcome, and one
#' variance parameter per outcome;
#' (7) the outcome-specific random intercepts u_ik: one variance parameter per outcome;
#' (8) the outcome-specific standard errors e_ijk: one variance parameter per outcome.
#' 
#' 
#' 
#' @param response character vector containing the names of the outcome-responses. ## Y
#' @param range vector indicating the range of the outcomes (that is
#' the minimum and maximum). ## Ak and Kk
#' @param var.time name of the time variable.   ## t
#' @param fixed_s a one-sided formula for specifying the fixed-effects defining 
#' the time-shift s.
#' 
#' @param subject name of the covariate representing the grouping structure.
#' @param data data frame containing all variables named in \code{response},
#'  \code{fixed_s}, \code{survival} and \code{subject}.
#' @param B optional specification for the initial values for the parameters.
#' Initial values should be entered in the order detailed in \code{details} section.
#' @param subset optional vector giving the subset of observations in
#' \code{data} to use. By default, all lines.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param posfix Optional vector giving the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all parameters are
#' estimated.
#' 
#' @param survival two-sided formula object. The left side of the formula corresponds 
#' to a Surv() object for right-censored (\code{Surv(EntryTime,Time,Indicator)})
#' and possibly left-truncated (\code{Surv(EntryTime,Time,Indicator)}).
#' Multiple causes of event can be considered 
#' in the Indicator (0 for censored, k for event of cause k). The right side of the 
#' formula specifies the covariates to include in the survival model.
#' Cause-specific covariate effect are specified with \code{cause()}.For example,  
#' Surv(Time,Indicator) ~ X1 + cause(X2) indicates a common effect of X1 and a cause-specific effect of X2. 
#' @param hazard optional family of hazard function assumed for the survival model. 
#' By default, "Weibull" specifies a Weibull baseline risk function. Other possibilities 
#' are "piecewise" for a piecewise constant risk function or "splines" for a cubic M-splines 
#' baseline risk function. For these two latter families, the number of nodes and the 
#' location of the nodes should be specified as well, separated by -. The number of 
#' nodes is entered first followed by -, then the location is specified with "equi", 
#' "quant" or "manual" for respectively equidistant nodes, nodes at quantiles of the 
#' times of event distribution or interior nodes entered manually in argument hazardnodes. 
#' It is followed by - and finally "piecewise" or "splines" indicates the family of 
#' baseline risk function considered. Examples include "5-equi-splines" for M-splines 
#' with 5 equidistant nodes, "6-quant-piecewise" for piecewise constant risk over 5 
#' intervals and nodes defined at the quantiles of the times of events distribution 
#' and "9-manual-splines" for M-splines risk function with 9 nodes, the vector of 7 
#' interior nodes being entered in the argument hazardnodes. In the presence of competing 
#' events, a vector of hazards should be provided such as hazard=c("Weibull","5-quant-splines") 
#' with 2 causes of event, the first one modelled by a Weibull baseline cause-specific 
#' risk function and the second one by splines.
#' @param hazardnodes optional vector containing interior nodes if splines or piecewise 
#' is specified for the baseline hazard function in hazard.
#' @param logscale optional boolean indicating whether an exponential
#' (logscale=TRUE) or a square (logscale=FALSE -by default) transformation is
#' used to ensure positivity of parameters in the baseline risk functions. See
#' details section
#' @param startWeibull optional numeric with Weibull hazard functions only.
#' Indicates the shift in the Weibull distribution.
#' @param hazardrange optional vector indicating the range of the survival times 
#' (that is the minimum and maximum). By default, the range is defined according 
#' to the minimum and maximum observed values of the survival times. The option 
#' should be used only for piecewise constant and Splines hazard functions.
#' 
#' @param convB optional threshold for the convergence criterion based on the
#' parameter stability. By default, convB=0.0001.
#' @param convL optional threshold for the convergence criterion based on the
#' log-likelihood stability. By default, convL=0.0001.
#' @param convG optional threshold for the convergence criterion based on the
#' derivatives. By default, convG=0.0001.
#' @param maxiter optional maximum number of iterations for the Marquardt
#' iterative algorithm. By default, maxiter=100.
#' @param nsim number of points used to plot the estimated link functions. By
#' default, nsim=100.    ## ?????
#' @param partialH optional vector giving the indices in vector B of parameters that
#' can be dropped from the Hessian matrix to define convergence criteria.
#' @param verbose logical indicating if information about computation should be
#' reported. Default to TRUE.
#' @param returndata logical indicating if data used for computation should be
#' returned. Default to FALSE, data are not returned.
#' @param methInteg character indicating the type of integration to compute the
#' log-likelihood. 'MCO' for ordinary Monte Carlo, 'MCA' for antithetic Monte Carlo,
#' 'QMC' for quasi Monte Carlo. Default to "QMC".
#' @param nMC integer, number of Monte Carlo simulations. Default to 1000.
#' @param nproc number of cores for parallel computation.
#' @param clustertype the type of cluster that should internally be created.
#' See \code{parallel::makeCluster} for possible values.
#' 
#' 
#' 
#' @return An object of class "LTSM" is returned containing some internal information
#' used in related functions. Users may investigate the following elements : 
#' \item{ns}{number of grouping units in the dataset} 
#' \item{loglik}{log-likelihood of the model} 
#' \item{best}{vector of parameter estimates in the same order as specified in 
#' \code{B} and detailed in section \code{details}}
#' \item{V}{vector containing the upper triangle matrix of variance-covariance
#' estimates of \code{best} with exception for variance-covariance parameters
#' of the random-effects for which \code{V} contains the variance-covariance
#' estimates of the Cholesky transformed parameters displayed in \code{cholesky}} 
#' \item{gconv}{vector of convergence criteria: 1. on the parameters, 2. on the 
#' likelihood, 3. on the derivatives} 
#' \item{conv}{status of convergence: =1 if the convergence criteria were satisfied, 
#' =2 if the maximum number of iterations was reached, =4 or 5 if a problem occured 
#' during optimisation} 
#' \item{call}{the matched call} 
#' \item{niter}{number of Marquardt iterations} 
#' \item{nevent}{number of occured event}
#' \item{pred}{table of individual predictions and residuals in the underlying
#' latent process scale; it includes marginal predictions (pred_m), marginal
#' residuals (resid_m), subject-specific predictions (pred_ss) and
#' subject-specific residuals (resid_ss) and finally the transformed
#' observations in the latent process scale (obs).}
#' \item{predRE_Y}{table containing individual predictions of the outcome-specific
#' random intercept}
#' \item{predSurv}{table containing the predicted baseline risk function and
#' the predicted cumulative baseline risk function }
#' \item{cholesky}{vector containing the estimates of the Cholesky transformed
#' parameters of the variance-covariance matrix of the random-effects}
#' \item{estimlink}{table containing the simulated values of each outcome and
#' the corresponding estimated link function} 
#' \item{AIC}{the Akaike's information criterion}
#' \item{BIC}{the Bayesian information criterion}
#' \item{CPUtime}{the runtime in seconds}
#' \item{data}{the original data set (if returndata is TRUE)}
#' @author Viviane Philipps, Tiphaine Saulnier and Cecile Proust-Lima
#' 
#' @references
#' Kühnel, Raket, Åström, Berger, Hansen, Krismer, ... & EMSA‐SG Natural History Study Investigators. (2022).
#' Disease Progression in Multiple System Atrophy—Novel Modeling Framework and Predictive Factors. Movement Disorders.
#' 
#' Philipps, Hejblum, Prague, Commenges, Proust-Lima (2021).
#' Robust and efficient optimization using a Marquardt-Levenberg algorithm with 
#' R package marqLevAlg, The R Journal 13:2.  
#' 
#' @examples
#' #### Examples
#' 
#' 
#' 
#' @export
#' 

################################################################################
################################################################################

LTSM <- function(response, range, var.time, fixed_s, 
    
                 subject,data,subset=NULL,na.action=1,B,posfix=NULL,
    
                 survival=NULL,hazard="Weibull",hazardrange=NULL,hazardnodes=NULL,logscale=FALSE,startWeibull=0, 
    
                 methInteg="QMC",nMC=1000,maxiter=100,convB=0.0001,convL=0.0001,convG=0.0001,partialH=NULL,
                 nsim=100,verbose=TRUE,returndata=FALSE,nproc=1, clustertype=NULL)
{
    ptm <- proc.time()
    
    cl <- match.call()
    
    if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
    if(nrow(data)==0) stop("Data should not be empty")
    if(!(na.action %in% c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")
    if(missing(subject)){ stop("The argument subject must be specified")}
    if(!(subject %in% colnames(data))) stop("Unable to find variable 'subject' in 'data'")
    if(!is.numeric(data[,subject])) stop("The argument subject must be numeric")
    
    nom.subject <- as.character(subject)
    
    if(missing(response)) stop("The argument response must be specified in any model")
    if(!all(response %in% colnames(data))) stop("Unable to find variable from 'response' in 'data'")
    if(missing(range)) stop("The argument range must be specified in any model")
    if(missing(var.time) | length(var.time)!=1)stop("The argument var.time is missing or is not of length 1")
    if(!(var.time %in% colnames(data))) stop("Unable to find variable 'var.time' in 'data'")
    if(missing(fixed_s)) stop("The argument fixed_s must be specified in any model")
    if(class(fixed_s)!="formula") stop("The argument fixed_s must be a formula")
    
    
    ## garder data tel quel pour le renvoyer
    if(returndata==TRUE)
    {
        datareturn <- data
    }
    else
    {
        datareturn <- NULL
    }
    
    
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
    
    ## argument subset
    form1 <- paste(c(nom.subject,nomsY,ttesLesVar),collapse="+")
    if(!isTRUE(all.equal(as.character(cl$subset),character(0))))
    {
        cc <- cl
        cc <- cc[c(1,which(names(cl)=="subset"))]
        cc[[1]] <- as.name("model.frame")
        cc$formula <- formula(paste("~",form1))
        cc$data <- data
        cc$na.action <- na.pass
        data <- eval(cc)
    }
    
    attributes(data)$terms <- NULL
    
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
    
    
    V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres
    
    ## parametres MC
    methInteg <- switch(methInteg,"MCO"=1,"MCA"=2,"QMC"=3)
    seqMC <- 0
    dimMC <- 0
    if(methInteg==3)
    {
        dimMC <- nvc.s + nvc.v + nvc.u   #+ ?
        if(dimMC>0) seqMC <- randtoolbox::sobol(n=nMC,dim=dimMC,normal=TRUE,scrambling=1)
    }
    
    
    ###valeurs initiales
    if(!(missing(B)))
    {
      if(!is.vector(B)) stop("B should be a vector")
      
      if (length(B)==NPM) b <- B
      else stop(paste("Vector B should be of length",NPM))
      
    }
    else ## B missing
    {
      b <- rep(0,NPM)
      
      if(nbevt>0){  #TS : si Weibull et prms au carre pr positivite -> valeurs par defaut = 1
        if(any(hazard!="Weibull")==FALSE & isFALSE(logscale)==TRUE){
          for(i in 1:nbevt){
            if(typrisq[i]==2){
              b[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- 1
              if(idtrunc==1 & any(Tentry==0)) #TS : si entree retardee et au moins un temps vaut 0
                b[sum(nprisq[1:i])-nprisq[i]+nprisq[i]] <- 1.25 #sinon pblm lors de recherche prms, risq instant tend vers infini qd w2 < 1 et t=0
            }
          }
        }
        if(any(hazard %in% c("splines","Splines")) & idtrunc==1 & any(Tentry==0)){
          for(i in 1:nbevt){
            if(typrisq[i]==3){
              b[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- 10^-7 #sinon pgrm crash
            }
          }
        }
      }
      
      #1 for variances
      b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s] <- 1
      b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + 1:nvc.v] <- rep(1,nmu.v)
      b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + 1:nvc.u] <- rep(1,nvc.u)
      b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + nvc.u + 1:nvc.err] <- rep(1,nvc.err)
      
    }
    
    
    ##------------------------------------------
    ##------nom au vecteur best
    ##--------------------------------------------
    
    nom.X0 <- colnames(X0)
    nom.X0[nom.X0=="(Intercept)"] <- "intercept"
    
    if(nbevt>0)
    {
      ##prm fct de risque
      if(isTRUE(logscale))
      {
        for(i in 1:nbevt)
        {
          nom1 <- rep(paste("event",i,sep=""),nprisq[i])
          if(typrisq[i]==2)
          {
            names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:2]," log(Weibull",1:2,")",sep="")
          }
          if(typrisq[i]==1)
          {
            names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(piecewise",1:(nz[i]-1),")",sep="")
          }
          if(typrisq[i]==3)
          {
            names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(splines",1:(nz[i]+2),")",sep="")
          }
        }
      }
      else
      {
        for(i in 1:nbevt)
        {
          nom1 <- rep(paste("event",i,sep=""),nprisq[i])
          if(typrisq[i]==2)
          {
            names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:2]," +/-sqrt(Weibull",1:2,")",sep="")
          }
          if(typrisq[i]==1)
          {
            names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(piecewise",1:(nz[i]-1),")",sep="")
          }
          if(typrisq[i]==3)
          {
            names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(splines",1:(nz[i]+2),")",sep="")
          }
        }
      }
      
      
      ##prm covariables survival
      nom1 <- NULL
      for(j in 1:nv)
      {
        if(idsurv[j]==1) #X
        {
          if(idtdv[j]==1)
          {
            nom1 <- c(nom1,paste("I(t>",nom.timedepvar,")",sep=""))
          }
          else
          {
            nom1 <- c(nom1,nom.X0[j])
          }
        }
        
        if(idsurv[j]==2) #cause(X)
        {
          if(idtdv[j]==1)
          {
            nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") event",1:nbevt,sep=""))
          }
          else
          {
            nom1 <- c(nom1,paste(nom.X0[j],paste("event",1:nbevt,sep="")))
          }
        }
        
      }
      
      if(nvarxevt>0) names(b)[nrisqtot+1:nvarxevt] <- nom1
      
      
      for(i in 1:nbevt)
        names(b)[nrisqtot+nvarxevt+(nbevt-1)*nea+1:nea] <- paste("event",i," asso",1:nea,sep="")
      
    }
    
    
    # partie long
    names(b)[nrisqtot+nvarxevt+nasso+1:nef.s] <- nom.X0[idg!=0]
    names(b)[nrisqtot+nvarxevt+nasso+nef.s+1] <- "var.s"
    names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+1:nef.B] <- paste("rate",1:ny,sep="")
    names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+1:nmu.v] <- paste("mu.v",1:ny,sep="")
    names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+1:nvc.v] <- paste("var.v",1:ny,sep="")
    names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+1:nvc.u] <- paste("var.u",1:ny,sep="")
    names(b)[nrisqtot+nvarxevt+nasso+nef.s+nvc.s+nef.B+nmu.v+nvc.v+nvc.u+1:nvc.err] <- paste("std.err",1:ny,sep="")
    
    namesb <- names(b)
    
    
    ## prm fixes
    fix <- rep(0,NPM)
    if(length(posfix))
    {
      if(any(!(posfix %in% 1:NPM))) stop("Indexes in posfix are not correct")
      
      fix[posfix] <- 1
    }
    if(length(posfix)==NPM) stop("No parameter to estimate")
    
    if(!all(partialH %in% 1:NPM)) stop(paste("partialH should contain indices between 1 and",NPM))
    
    # RE : var -> std
    b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s] <- sqrt(b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s])
    b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + 1:nvc.v] <- sqrt(b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + 1:nvc.v])
    b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + 1:nvc.u] <- sqrt(b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + 1:nvc.u])
    #b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + nvc.u + 1:nvc.err] <- sqrt(b[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + nvc.u + 1:nvc.err])
    
    
    # fixed prms
    nfix <- sum(fix)
    bfix <- 0
    if(nfix>0)
    {
      bfix <- b[which(fix==1)]
      b <- b[which(fix==0)]
      NPM <- NPM-nfix
    }
    
    
    #browser()
    ###estimation
    
    conv3 <- c(convB, convL, convG)
    
    
    ###############
    ###   MLA   ###
    ###############
    
    if(maxiter==0)
    {
      vrais <- loglik(b,Y0,X0,idg,
                        ny,ns,nv,
                        nobs,nmes,NPM,nfix,bfix,
                        zitr,fix,
                        methInteg,nMC,dimMC,seqMC)

      out <- list(conv=2, V=rep(NA, length(b)), best=b, predRE_Y=NA,
                  Yobs=NA, resid_m=NA, resid_ss=NA, risqcum_est=NA, risq_est=NA,
                  marker=NA, transfY=NA, gconv=rep(NA,3), niter=0, loglik=vrais)
    }
    else
    {
      res <- mla(b=b, m=length(b), fn=loglik,
                 clustertype=clustertype,.packages=NULL,
                 epsa=convB,epsb=convL,epsd=convG,
                 digits=8,print.info=verbose,blinding=FALSE,
                 multipleTry=25,file="",partialH=partialH,
                 nproc=nproc,maxiter=maxiter,minimize=FALSE,
                 Y0=Y0,X0=X0,idg0=idg,
                 ny0=ny,ns0=ns,nv0=nv,
                 nobs0=nobs,nmes0=nmes,npm0=NPM,nfix0=nfix,bfix0=bfix,
                 zitr0=zitr,fix0=fix,
                 methInteg0=methInteg,nMC0=nMC,dimMC0=dimMC,seqMC0=seqMC)

      out <- list(conv=res$istop, V=res$v, best=res$b, predRE_Y=NA, Yobs=NA,
                  resid_m=NA, resid_ss=NA, risqcum_est=NA, risq_est=NA, marker=NA,
                  transfY=NA, gconv=c(res$ca, res$cb, res$rdm), niter=res$ni,
                  loglik=res$fn.value)
    }
    
    
    ## creer best a partir de b et bfix
    best <- rep(NA,length(fix))
    best[which(fix==0)] <- out$best
    best[which(fix==1)] <- bfix
    out$best <- best
    NPM <- NPM+nfix


    ## mettre NA pour les variances et covariances non calculees et 0 pr les prm fixes
    if(length(posfix))
    {
      mr <- NPM-length(posfix)
      Vr <- matrix(0,mr,mr)
      Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
      Vr <- t(Vr)
      Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
      V <- matrix(0,NPM,NPM)
      V[setdiff(1:NPM,posfix),setdiff(1:NPM,posfix)] <- Vr
      V <- V[upper.tri(V,diag=TRUE)]
    }
    else
    {
      V <- out$V
    }


    ## RE : std -> var
    out$best[nrisqtot + nvarxevt + nasso + nef.s + nvc.s] <- out$best[nrisqtot + nvarxevt + nasso + nef.s + nvc.s]**2
    out$best[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + 1:nvc.v] <- out$best[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + 1:nvc.v]**2
    out$best[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + 1:nvc.u] <- out$best[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + 1:nvc.u]**2
    #out$best[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + nvc.u + 1:nvc.err] <- out$best[nrisqtot + nvarxevt + nasso + nef.s + nvc.s + nef.B + nmu.v + nvc.v + nvc.u + 1:nvc.err]**2  # std


    ##pred
    pred_m <- out$Yobs-out$resid_m
    pred_ss <- out$Yobs - out$resid_ss
    pred <- data.frame(IND,outcome,pred_m,out$resid_m,pred_ss,out$resid_ss,out$Yobs)

    colnames(pred)<-c(nom.subject,"Yname","pred_m","resid_m","pred_ss","resid_ss","obs")
    rownames(pred) <- NULL

    ## risques
    if(nbevt>0)
    {
      risqcum_est <- matrix(out$risqcum_est,nrow=nsim,ncol=nbevt)
      risq_est <- matrix(out$risq_est,nrow=nsim,ncol=nbevt)
      predSurv <- cbind(time,risq_est,risqcum_est)

      temp <- paste("event",1:nbevt,".RiskFct",sep="")
      temp1 <- paste("event",1:nbevt,".CumRiskFct",sep="")
      colnames(predSurv) <- c("time",temp,temp1)
      rownames(predSurv) <- 1:nsim
    }
    else
    {
      predSurv <- NA
    }


    N <- rep(NA,14)
    N[1] <- 0 #nprob
    N[2] <- nrisqtot
    N[3] <- nvarxevt
    N[4] <- nasso
    N[5] <- nef.s
    N[6] <- nvc.s
    N[7] <- nef.B
    N[8] <- nmu.v
    N[9] <- nvc.v
    N[10] <- nvc.u
    N[11] <- nvc.err
    N[12] <- ny
    N[13] <- nobs
    N[14] <- nbevt

    nevent <- rep(0,nbevt)
    for(ke in 1:nbevt)
    {
      nevent[ke] <- length(which(devt==ke))
    }

    Nprm <- c(nprisq)


    nom.X0[nom.X0=="(Intercept)"] <- "Intercept"



    ## noms des variables
    Names <- list(Xnames=nom.X0[-1],Ynames=nomsY,
                  ID=nom.subject,Tnames=noms.surv,
                  Xvar=setdiff(ttesLesVar,noms.surv))

    form <- list(fixed_s=fixed_s,
                 form.commun=form.commun, form.cause=form.cause)

    cost <- proc.time()-ptm

    res <-list(ns=ns,idg=idg,idsurv=idsurv,
               loglik=out$loglik,best=out$best,V=V,gconv=out$gconv,conv=out$conv,
               call=cl,niter=out$niter,N=N,nevent=nevent,Nprm=Nprm,
               pred=pred,Names=Names,form=form,
               logspecif=logspecif,predSurv=predSurv,typrisq=typrisq,hazardnodes=zi,nz=nz,
               linknodes=zitr,
               na.action=nayk,AIC=2*(length(out$best)-length(posfix)-out$loglik),BIC=(length(out$best)-length(posfix))*log(ns)-2*out$loglik,data=datareturn,
               #wRandom=wRandom,b0Random=b0Random,
               posfix=posfix,CPUtime=cost[3])

    names(res$best) <- namesb
    class(res) <-c("LTSM")
    
    
    return(res)
    
    
}

loglik <- function(b0,Y0,X0,idg0,
                       ny0,ns0,nv0,
                       nobs0,nmes0,npm0,nfix0,bfix0,
                       zitr0,fix0,
                       methInteg0,nMC0,dimMC0,seqMC0)
{
    res <- 0
    .Fortran(C_loglik,as.double(Y0),as.double(X0),as.integer(idg0),as.integer(ny0),as.integer(ns0),as.integer(nv0),as.integer(nobs0),as.integer(nmes0),as.integer(npm0),as.double(b0),as.integer(nfix0),as.double(bfix0),as.double(zitr0),as.integer(fix0),as.integer(methInteg0),as.integer(nMC0),as.integer(dimMC0),as.double(seqMC0),loglik_res=as.double(res))$loglik_res
}
