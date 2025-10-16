#' Estimation of latent time shift model for multivariate longitudinal 
#' outcomes (and time-to-event data) with a time-shifting process.
#' 
#' This function fits extended joint disease progression models with shared random effects.
#' The aim of the joint disease progression model is to describe the markers’ trajectories
#' over the latent disease time, reconstituted by relocating patients along the temporal
#' axis according to patient’s disease advancement and marker observations, while
#' simultaneously accounting for the association with the risk of event occurence during
#' follow-up. The model consists of three submodels : one for the temporal recalibration,
#' one for the longitudinal markers, and one for the time-to-event process.
#' The individual temporal recalibration results from observed time regressed on covariates,
#' and an individual deviation.
#' The longitudinal submodel handles multiple continuous longitudinal outcomes
#' (Gaussian or non-Gaussian, curvilinear) in a mixed effects framework.
#' The survival submodel handles right-censored time-to-events with competing risks.
#' The association between the longitudinal and the survival data is captured by including 
#' the random effect from the mixed model or the predicted current level 
#' of the underlying process as a linear predictor in the proportional hazard survival model.
#' Parameters of the submodels are estimated simultaneously using a maximum likelihood method,
#' through a Marquardt-Levenberg algorithm.
#' 
#' 
#' @param fixed a list containing two-sided linear formula object for specifying the
#' fixed-effects in the marker-specific submodel. The response outcomes are separated by 
#' \code{+} on the left of \code{~} and the covariates are separated by \code{+} on 
#' the right of the \code{~}. For data-driven trajectory, polynomial time functions 
#' should be specified with \code{FP()} function, with powers specified between brackets. 
#' For sigmoid trajectory, no time function should be specified.
#' @param random a list containing one-sided formula for the random-effects in the
#' marker-specific submodel. Covariates with a random-effect are separated
#' by \code{+}. For sigmoid trajectory, NULL should be specified.
#' @param range a list containing vector indicating the range of the outcomes (that is
#' the minimum and maximum).
#' @param link optional vector of families of parameterized link functions defining
#' the measurement models (one by outcome). Option "id" (default) indicates no link function.
#' Other possibilities include "linear" specifies a linear link function, and "Splines" for 
#' approximating the link function by I-splines. 
#' For splines case, the number of nodes and the nodes location should be also 
#' specified. The number of nodes is first entered followed by \code{-}, then the location 
#' is specified with "equi", "quant" or "manual" for respectively equidistant
#' nodes, nodes at quantiles of the marker distribution or interior nodes
#' entered manually in argument \code{intnodes}. It is followed by \code{-splines}.
#' For example, "7-equi-splines" means
#' I-splines with 7 equidistant nodes, "6-quant-splines" means I-splines with 6
#' nodes located at the quantiles of the marker distribution and
#' "9-manual-splines" means I-splines with 9 nodes, the vector of 7 interior
#' nodes being entered in the argument \code{intnodes}.
#' @param intnodes optional vector of interior nodes. This argument is only
#' required for a I-splines link function with nodes entered manually.
#' @param model vector defining the model basis (one per outcome). Option "logistic" indicates
#' a sigmoid trajectory model, and "linear" indicates a data-driven trajectory model.
#' @param timevar name of the variable representing the measurement times.
#' @param fixed.s one-sided formula for covariates in the
#' temporal recalibration submodel. Covariates are separated
#' by \code{+}.
#' @param subject name of the covariate representing the grouping structure.
#' @param data data frame containing all variables named in other arguments.
#' @param subset optional vector giving the subset of observations in
#' \code{data} to use. By default, all lines.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param B optional specification for the initial values for the parameters.
#' @param posfix Optional vector giving the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all parameters are
#' estimated.
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
#' @param hazardrange optional vector indicating the range of the survival times 
#' (that is the minimum and maximum). By default, the range is defined according 
#' to the minimum and maximum observed values of the survival times. The option 
#' should be used only for piecewise constant and Splines hazard functions.
#' @param hazardnodes optional vector containing interior nodes if splines or piecewise 
#' is specified for the baseline hazard function in hazard.
#' @param logscale optional boolean indicating whether an exponential
#' (logscale=TRUE) or a square (logscale=FALSE -by default) transformation is
#' used to ensure positivity of parameters in the baseline risk functions. See
#' details section.
#' @param startWeibull optional numeric with Weibull hazard functions only.
#' Indicates the shift in the Weibull distribution.
#' @param surv.shift optional one-sided formula specifying the time function in interaction to the 
#' individual shift in the survival submodel. Default is NULL.
#' @param sharedtype vector indicating for each outcome and each event the shared random function type. 
#' \code{'RE'} indicates an association through the random effects included in the marker-specific submodel.
#' \code{'CV'} defines a association through the predicted current value of the marker, in the latent process scale.
#' @param nGK Number of integration points from Gauss-Kronrod procedure (for association through marker current value).
#' Possibilities are 7 or 15 (default).
#' @param methInteg character indicating the type of integration to compute the
#' log-likelihood. 'MCO' for ordinary Monte Carlo, 'MCA' for antithetic Monte Carlo,
#' 'QMC' for quasi Monte Carlo. Default to "QMC".
#' @param nMC integer, number of Monte Carlo simulations. Default to 1000.
#' @param maxiter optional maximum number of iterations for the Marquardt
#' iterative algorithm. By default, maxiter=100.
#' @param convB optional threshold for the convergence criterion based on the
#' parameter stability. By default, convB=0.0001.
#' @param convL optional threshold for the convergence criterion based on the
#' log-likelihood stability. By default, convL=0.0001.
#' @param convG optional threshold for the convergence criterion based on the
#' derivatives. By default, convG=0.0001.
#' @param partialH optional vector giving the indices in vector B of parameters that
#' can be dropped from the Hessian matrix to define convergence criteria.
#' @param nsim number of points used to plot the estimated link functions. By
#' default, nsim=100.
#' @param verbose logical indicating if information about computation should be
#' reported. Default to TRUE.
#' @param returndata logical indicating if data used for computation should be
#' returned. Default to FALSE, data are not returned.
#' @param nproc number of cores for parallel computation.
#' @param clustertype the type of cluster that should internally be created.
#' See \code{parallel::makeCluster} for possible values.
#' 
#' 
#' @return An object of class "jointLTSM".
#' 
#' @author Tiphaine Saulnier, Cecile Proust-Lima
#' 
#' @references
#' Saulnier T., Fabbri M., Pavy-Le Traon A., Le Goff M., Helmer C., Péran P., Meissner W. G., 
#' Rascol O., Foubert-Samier A., Proust-Lima C. (2023). Disease Progression in Multiple System 
#' Atrophy : The Value of Clinical Cohorts with Long Follow-Up. Movement Disorders, 38(8), 1567-1569.
#' 
#' @export
#' 

################################################################################
################################################################################

jointLTSM <- function(fixed, random, range,
                      link, intnodes=NULL,
                      model,
                      timevar,
                      
                      fixed.s,
                      
                      subject,data,subset=NULL,na.action=1,B,posfix=NULL,
                      
                      survival=NULL,hazard="Weibull",hazardrange=NULL,hazardnodes=NULL,logscale=FALSE,startWeibull=0, 
                      surv.shift=NULL,
                      sharedtype="RE", nGK=15,
                      
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
    
    if(missing(fixed)) stop("The argument fixed must be specified in any model")
    if(!inherits(fixed,"list")) stop("The argument fixed must be a list of 2-sided formula(s)")
    nlp <- length(fixed) # nb of latent processes
    if(missing(random)) stop("The argument random must be specified in any model")
    if(!inherits(random,"list")) stop("The argument random must be a list of 1-sided formula(s)")
    if(length(random)!=nlp) stop(paste("The argument random should be a list of length ",nlp,sep=""))
    if(missing(range)) range <- rep(list(NULL),nlp)
    if(!inherits(range,"list")) stop("The argument range must be a list")
    if(length(range)!=nlp) stop(paste("The argument range should be a list of length ",nlp,sep=""))
    if(missing(link)) stop("The argument link must be specified in any model")
    if(length(link)!=nlp) stop(paste("The argument link should be a vector of length ",nlp,sep=""))
    if(all(link %in% c("id","linear")) & !is.null(unlist(intnodes))) stop("Intnodes should only be specified with splines links")
    if(missing(model)) stop("The argument model must be specified in any model")
    if(length(model)!=nlp) stop(paste("The argument model should be a vector of length ",nlp,sep=""))
    if(!any(model %in% c("logistic","linear"))) stop("The argument model must only contain 'logistic', or 'linear'")
    if(missing(timevar) | length(timevar)!=1) stop("The argument timevar is missing or is not of length 1")
    if(!(timevar %in% colnames(data))) stop("Unable to find variable 'timevar' in 'data'")
    
    if(missing(fixed.s)) stop("The argument fixed.s must be specified")
    if(class(fixed.s)!="formula") stop("The argument fixed.s must be a formula")
    
    
    ## garder data tel quel pour le renvoyer
    if(returndata==TRUE)
    {
        datareturn <- data
    }
    else
    {
        datareturn <- NULL
    }
    
    
    ##### liste par latent process
    ttesLesVar <- NULL
    nomsY <- NULL
    ny <- 0
    for(l in 1:nlp){
      
      ##pour acces aux attributs des formules
      # fixed
      if(!inherits(fixed[[l]],"formula")) stop("The argument fixed must contain 2-sided formula(s)")
      afixed <- terms(fixed[[l]], specials=c("factor"))
      fixed2 <- formula(paste("~",fixed[[l]][3],sep=""))
      ttesLesVar <- c(ttesLesVar,colnames(get_all_vars(terms(fixed2),data=data[1,])))
      
      # random
      if(!is.null(random[[l]])){
        if(!inherits(random[[l]],"formula") & !is.null(random[[l]])) stop("The argument random must be NULL or contain 1-sided formula(s)")
        arandom <- terms(random[[l]], specials=c("factor"))
        ttesLesVar <- c(ttesLesVar,colnames(get_all_vars(terms(random[[l]]),data=data[1,])))
      }
      
      ##verifier si toutes les variables sont dans data
      variables <- c(attr(afixed,"variables"))
      if(!is.null(random[[l]])) variables <- c(variables,attr(arandom,"variables"))
      variables <- unlist(lapply(variables,all.vars))  
      if(!all(variables %in% colnames(data))) stop(paste("Data should contain the variables",paste(unique(variables),collapse=" ")))
      
      ##liste des outcomes
      nomsYl <- as.character(attr(afixed,"variables")[2])
      nomsY <- c(nomsY,nomsYl)
      # nomsY <- strsplit(nomsY,split=" + ",fixed=TRUE)
      # nomsY <- as.vector(nomsY[[1]])
      nyl <- length(nomsYl)
      ny <- ny + nyl
      if(nyl > 1) stop("dimension process not coded yet")
      
    }
    
    
    ## formule du time shift
    ttesLesVar <- c(ttesLesVar,colnames(get_all_vars(terms(fixed.s),data=data[1,])))
    
    ##verifier si toutes les varialbes sont dans data
    ttesLesVar <- unique(ttesLesVar)
    #if(any(ttesLesVar=="FP")) stop("For code purpose, no variable can be named 'FP'.")
    if(!all(ttesLesVar %in% colnames(data))) stop(paste("Data should contain the variables",paste(ttesLesVar,collapse=" ")))
    
    
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
      typrisq <- 0
      noms.surv <- NULL
      form.cause <- ~-1
      surv.shift <- NULL
      survival <- ~-1
      nz <- 0
      zi <- 0
      minT <- 0
      maxT <- 0
    }
    else
    {
      
      ## objet Surv
      surv <- cl$survival[[2]] #surv <- survival[[2]]
      
      if(length(surv)==3) #censure droite sans troncature gauche
      {
        idtrunc <- 0
        
        Tevent <- getElement(object=data,name=as.character(surv[2]))
        Event <- getElement(object=data,name=as.character(surv[3]))
        Tentry <- rep(0,length(Tevent)) #si pas de troncature, Tentry=0
        
        noms.surv <-  c(as.character(surv[2]),as.character(surv[3]))
        
        surv <- do.call("Surv",list(time=Tevent,event=Event,type="mstate"))
      }
      
      if(length(surv)==4) #censure droite et troncature
      {
        idtrunc <- 1
        
        Tentry <- getElement(object=data,name=as.character(surv[2]))
        Tevent <- getElement(object=data,name=as.character(surv[3]))
        Event <- getElement(object=data,name=as.character(surv[4]))
        
        noms.surv <-  c(as.character(surv[2]),as.character(surv[3]),as.character(surv[4]))
        
        surv <- do.call("Surv",list(time=Tentry,time2=Tevent,event=Event,type="mstate"))
      }
      
      ## nombre d'evenement concurrents
      nbevt <- length(attr(surv,"states"))
      if(nbevt<1) stop("No observed event in the data")
      
      
      ## pour la formule pour survival
      form.surv <- cl$survival[3] #form.surv <- survival[3]
      
      noms.form.surv <- all.vars(attr(terms(formula(paste("~",form.surv))),"variables"))
      if(length(noms.form.surv)==0)
      {
        form.cause <- ~-1
        asurv <- terms(~-1)
      }
      else
      {
        ##creer la formula
        # par defaut les effets sont specifiques aux events
        form1 <- formula(paste("~",form.surv))
        asurv <- terms(form1)
        ind.cause <- 1:length(labels(asurv))
        if(length(ind.cause))
        {
          form.cause <- paste(labels(asurv)[ind.cause],collapse="+")
          form.cause <- formula(paste("~",form.cause))
        }
        else
        {
          form.cause <- ~-1
        }
      }
      
      ##verifier si toutes les variables sont dans data
      varSurv <- unique(all.vars(terms(survival)))
      if(!all(varSurv %in% colnames(data))) stop(paste("Data should contain the variables",paste(varSurv,collapse=" ")))
      varSurv <- varSurv[!(varSurv %in% noms.surv)]
      ttesLesVar <- unique(c(ttesLesVar,varSurv))
      
      
      
      ## pour la formule pour surv.shift
      if(is.null(surv.shift))  surv.shift <- ~ 1
      if(length(surv.shift)!=2) stop("Argument surv.shift should be a one-sided formula")
      
      varSurvShift <- all.vars(terms(surv.shift))
      if(length(varSurvShift)) stop("Argument surv.shift can not contain covariates, only functions of time via FP")
      
    }
    
    ttesLesVar <- c(timevar,ttesLesVar)
    
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
    
    
    # // note : inversement paragraphes creation X0 et transfoY
    # car structure de X0 dependant de idmodel
    
    ## transfo outcomes
    
    # argument link
    if (length(link)!=1 & length(link)!=ny) stop("One link per outcome should be specified")
    if(any(link %in% c("splines","Splines"))) # par defaut : splines = 5-quant-splines
    {
      link[which(link %in% c("splines","Splines"))] <- "5-quant-splines"
    }
    if(length(link)==1 & ny>1)
    {
      link <- rep(link, ny)
    }
    
    idlink <- rep(2,ny) # splines
    idlink[which(link == "id")] <- 0
    idlink[which(link == "linear")] <- 1
    
    spl <- strsplit(link[which(idlink==2)],"-")
    if(any(sapply(spl,length)!=3)) stop("Invalid argument 'link'")
    
    nySPL <- length(spl)
    
    ##remplir range si pas specifie
    if(!is.null(range) & length(range)!=ny) stop("Length of list range is not correct.")
    rg <- NULL
    for(k in 1:ny){
      rg_obs <- range(get(dataY[k])[,nomsY[k]])
      rg_thq <- rg_obs
      if(!is.null(range[[k]])){
        rg_thq <- range[[k]][1:2]
        if(length(rg_thq)!=2) stop("In the argument range, a vector of 2 values (min and max) should be mentionned.")
        if(rg_obs[1]<rg_thq[1] | rg_obs[2]>rg_thq[2]) stop(paste("For marker ",nomsY[k],", the range specified do not cover the entire range of the data"),sep="")
      }
      rg <- rbind(rg,rg_thq)
    }
    row.names(rg) <- NULL
    
    
    # range et noeuds internes
    nbzitr <- rep(2,ny) #nbzitr = nb de noeuds si splines, 2 sinon (pr min/max)
    nbnodes <- NULL  #que pour les splines
    spltype <- NULL
    if(nySPL>0)
    {
      for (i in 1:nySPL)
      {
        nbnodes <- c(nbnodes, spl[[i]][1])
        spltype <- c(spltype, spl[[i]][2])
        if(spl[[i]][3] != "splines") stop("Invalid argument link")
      }
    }
    nbnodes <- as.numeric(nbnodes)
    nbzitr[which(idlink==2)] <- nbnodes
    
    ##test splines
    if(!(all(spltype %in% c("equi","quant","manual")))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")
    
    ##tester longueur de intnodes
    if(!is.null(unlist(intnodes)))
    {
      if(length(unlist(intnodes)) != sum(nbnodes[which(spltype=="manual")]-2)) stop(paste("Argument intnodes should be list of total length",sum(nbnodes[which(spltype=="manual")]-2)))
    }
    
    ##intnodes2 : contient tous les noeuds interieurs (pas seulement ceux de manual)
    intnodes2 <- rep(NA,sum(nbnodes-2))
    nb <- 0
    nbspl <- 0
    for (k in 1:ny)
    {
      if (idlink[k]!=2) next
      else
      {
        nbspl <- nbspl+1
        
        if(spltype[nbspl]=="manual")
        {
          nodes <- intnodes[[k]]
          if(!length(nodes)) stop("The length of intnodes is not correct")
          intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <-  nodes
          nb <- nb+nbnodes[nbspl]-2
          
          if(any(nodes <= rg[k,1]) | any(nodes >= rg[k,2])) stop("Interior nodes must be in the range of the outcome")
        }
        
        if(spltype[nbspl]=="equi")
        {
          nodes <- seq(range[2*(nbspl-1)+1], range[2*nbspl], length.out=nbnodes[nbspl])
          nodes <- nodes[-nbnodes[nbspl]]
          nodes <- nodes[-1]
          intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- nodes
        }
        
        if(spltype[nbspl]=="quant")
        {
          nodes <- quantile(get(dataY[k])[,nomsY[k]], probs=seq(0,1,length.out=nbnodes[nbspl]))
          if(length(unique(nodes)) != length(nodes)) stop(paste("Some nodes are equal for link number",k,"; Please try to reduce the number of nodes or use manual location."))
          nodes <- nodes[-nbnodes[nbspl]]
          nodes <- nodes[-1]
          intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- as.vector(nodes)
        }
      }
    }
    
    if(nb != length(unlist(intnodes))) stop(paste("The vector intnodes should be of length",nb))

    
    ##remplir zitr
    m <- 0
    if(nySPL>0) m <- max(nbnodes)
    zitr <- matrix(0,max(m,2),ny)
    nb12 <- 0
    nbspl <- 0
    for (k in 1:ny)
    {
      
      if((idlink[k]==0) | (idlink[k]==1)){
        nb12 <- nb12 + 1
        zitr[1:2,k] <- rg[k,]
      } 
      
      if(idlink[k]==2)
      {
        nb12 <- nb12+1
        nbspl <- nbspl+1
        zitr[2:(nbzitr[k]-1),k] <- intnodes2[ifelse(nbspl==1,0,1)*sum(nbnodes[1:(nbspl-1)]-2) + 1:(nbnodes[nbspl]-2)]
        zitr[1,k] <- rg[k,1]
        zitr[nbnodes[nbspl],k]  <- rg[k,2]
        
        ##verifier s'il y a des obs entre les noeuds
        hcounts <- hist(get(dataY[k])[,nomsY[k]],breaks=zitr[1:nbnodes[nbspl],k],plot=FALSE,include.lowest=TRUE,right=TRUE)$counts
        if(any(hcounts==0)) stop(paste("Link function number",k,"can not be estimated. Please try other nodes such that there are observations in each interval."))
      }
    }
    
    ##uniqueY0 et indiceY0
    uniqueY0 <- NULL
    indiceY0 <- NULL
    nvalSPLORD <- rep(0,ny)
    #modalites <- vector("list",ny)
    nb <- 0
    for (k in 1:ny)
    {
      
      if(idlink[k]!=2)
      {
        indiceY0 <- c(indiceY0, rep(0,length(get(dataY[k])[,nomsY[k]])))
        next
      }
      
      yk <- get(dataY[k])[,nomsY[k]]
      uniqueTemp <- sort(unique(yk))
      permut <- order(order(yk))  # sort(y)[order(order(y))] = y
      
      if(length(as.vector(table(yk)))==length(uniqueTemp))
      {
        indice <- rep(1:length(uniqueTemp), as.vector(table(yk)))
        if(idlink[k]==2)
        {
          indiceTemp <- nb + indice[permut]
        }
        else
        {
          indiceTemp <- indice[permut]
        }
        
        nb <- nb + length(uniqueTemp)
        
        uniqueY0 <- c(uniqueY0, uniqueTemp)
        indiceY0 <- c(indiceY0, indiceTemp)
        nvalSPLORD[k] <- length(uniqueTemp)
      }
      else
      {
        uniqueY0 <- c(uniqueY0, yk)
        indiceY0 <- c(indiceY0, ifelse(idlink[k]==1,nb,0)+c(1:length(yk)))
        nb <- nb + length(yk)
        nvalSPLORD[k] <- length(yk)
      }
      
    }
    
    
    ## outcome models
    # underlying process or trajectory
    
    # argument model
    if (length(model)!=1 & length(model)!=ny) stop("One model per outcome should be specified")
    if(length(model)==1 & ny>1)
    {
      model <- rep(model, ny)
    }
    
    idmodel <- rep(1,ny) # linear
    idmodel[which(model == "logistic")] <- 0

    
    # verification et reconstruction des elements des formules 
    # (partie longitudinale)
    fixed.Y <- NULL
    FP.power.fx <- NULL
    fixed.Y.interac <- NULL
    random.Y <- NULL
    FP.power.rd <- NULL
    for(k in 1:ny){
      
      fixed.Y.k <- fixed[[k]]
      random.Y.k <- random[[k]]
      
      if(timevar %in% all.vars(fixed.Y.k)) stop("In fixed, variable timevar should not appear in the formula.")
      if(timevar %in% all.vars(random.Y.k)) stop("In random, variable timevar should not appear in the formula.")
      if(!is.null(random.Y.k)){if(random.Y.k == "~-1") stop("Argument random should at least specify one.")}
      
      
      if(idmodel[k]==0){ # logistic
        
        # fixed
        if(gsub("contrast","",fixed.Y.k)[3]!="1") fixed.Y <- c(fixed.Y, list(formula(paste("~",gsub("contrast","",fixed.Y.k)[3]))))
        else fixed.Y <- c(fixed.Y, list(formula("~1")))
        FP.power.fx <- c(FP.power.fx,list(NULL))
        fixed.Y.interac <- c(fixed.Y.interac,0)
        
        # random
        if(!is.null(random.Y.k)){if(random.Y.k!="~1") stop("In random and for a logistic model, only one random intercept is implemented, the random structure cannot be specified, argument should be NULL or '~1'.")}
        #random.Y <- c(random.Y, list(formula(~1)))
        random.Y <- c(random.Y, 1)
        FP.power.rd <- c(FP.power.rd,list(NULL))
        
      }
      
      if(idmodel[k]==1){ # linear
        
        # fixed
          # time fct
        if(!("FP" %in% all.names(fixed.Y.k,unique=T))) stop("In fixed and for a linear model, time function should be specified via FP(), with powers between parentheses.")
        FP.power.fx <- c(FP.power.fx,list(as.numeric(unlist(strsplit(gsub('.*FP\\((.*)\\).*', '\\1', attr(terms(fixed.Y.k)[1],"term.labels")),split=",")))))
          # covariates
        var <- all.vars(fixed.Y.k)[-1]
        if(idlink[k]==0) var <- c(1,var) # intercept fixe
        if(length(var)) fixed.Y <- c(fixed.Y, list(formula(paste("~",paste(var,collapse="+")))))
        else fixed.Y <- c(fixed.Y, list(NULL))
          # interaction
        if(":" %in% all.names(fixed.Y.k,unique=T)){
          fixed.Y.interac <- c(fixed.Y.interac,1)
          print("By default, for a linear model, if an interaction term appears in fixed, all covariates will be in interaction with the time function. You can choose to remove an interaction effect by fixing the associated parameter to 0 using B and posfix arguments.")
          # note : si interaction entre deux covariables : creer la variable interaction au prealable
        }
        else fixed.Y.interac <- c(fixed.Y.interac,0)
        
        # random
          # time fct
        if(!("FP" %in% all.names(random.Y.k,unique=T)) & random.Y.k != "~1") stop("In random and for a linear model, time function should be specified via FP(), with powers between parentheses.")
        if(random.Y.k != "~1"){
          FP.power.rd <- c(FP.power.rd,list(as.numeric(unlist(strsplit(gsub('.*FP\\((.*)\\).*', '\\1', attr(terms(random.Y.k),"term.labels")),split=",")))))
        }
          # covariates
        if(length(all.vars(random.Y.k)[-1])) stop("For a linear model, no covariate possible in random.")
          # interaction
        if(":" %in% all.names(random.Y.k,unique=T)) stop("For a linear model, no interaction possible between time function and covariates in random.")
          # intercept
        if(attr(terms(random.Y.k),"intercept")==1)
          random.Y <- c(random.Y, 1)
        else random.Y <- c(random.Y, 0)

      }
      
    }
    
    
    # (partie shift)
    
      # long
    if(timevar %in% all.vars(fixed.s)) stop("In fixed.s, variable timevar can not appear in the formula.")
    
    # verif non time dependant covariate(s)
    if(length(all.vars(fixed.s))){
      if(nrow(unique(na.omit(data0[,c(subject,all.vars(fixed.s))]))) != length(unique(IND)))
        stop("In fixed.s, covariates should be time-independent.")
    }
    
      # surv
    surv.s <- NULL
    FP.power.s <- NULL
    
    if(nbevt > 0){
      
      # time fct
      if(!("FP" %in% all.names(surv.shift,unique=T)) & surv.shift != "~1") stop("In surv.shift, time function should be specified via FP(), with powers between parentheses.")
      if(surv.shift != "~1"){
        FP.power.s <- c(FP.power.s,list(as.numeric(unlist(strsplit(gsub('.*FP\\((.*)\\).*', '\\1', attr(terms(surv.shift),"term.labels")),split=",")))))
      }
      if(attr(terms(surv.shift),"intercept")==1)
        surv.s <- 1
      else surv.s <- 0
    
    }  
    
    
    ##creation de X0 (ttes les var)
    
    # vecteur des temps
    Xtime <- model.matrix(formula(paste("~",timevar)), data=data0)[,2:1] # time and intercept
    
    # matrices varexp pour les outcomes Y
    Xfixed.Y <- NULL
    z.fixed.Y <- NULL
    Xrandom.Y <- NULL
    z.random.Y <- NULL
    for(k in 1:ny){
      
      fixed.Y.k <- fixed.Y[[k]]
      random.Y.k <- random.Y[k]
      
      ## fixed effects
      if(!is.null(fixed.Y.k)){
        Xfixed.Y.k <- model.matrix(fixed.Y.k, data=data0)
        if(idlink[k]!=0 & idmodel[k]==1) Xfixed.Y.k <- Xfixed.Y.k[,-1,drop=F] # contrainte identifiabilite : no fixed intercept
        Xfixed.Y <- cbind(Xfixed.Y,Xfixed.Y.k)
        z.fixed.Y <- c(z.fixed.Y,list(colnames(Xfixed.Y.k)))
      }
      else z.fixed.Y <- c(z.fixed.Y,list(NULL))
      
      ## random effects 
      z.random.Y <- c(z.random.Y,random.Y.k)
      
    }
    Xfixed.Y <- Xfixed.Y[,unique(colnames(Xfixed.Y)),drop=F]
    
    # matrices varexp pour le time shift s
    
    Xfixed.s <- model.matrix(fixed.s, data=data0)
    
    # verif non time dependant covariate(s)
    if(length(all.vars(fixed.s))) Xfixed.s <- Xfixed.s[,-1,drop=F] # sans intercept si covariables
    z.fixed.s <- colnames(Xfixed.s)
    
    
    # matrices de survie
    
    # effets par cause
    Xsurvcause <- model.matrix(form.cause,data=data0)
    
    if(form.cause != ~-1)
    {
      z.survcause <- strsplit(colnames(Xsurvcause),split=":",fixed=TRUE)
      z.survcause <- lapply(z.survcause,sort)
    }
    else
    {
      z.survcause <- list() 
    }
    
    # fct du temps en interaction avec shift
    
    Xsurv.s <- NULL 
    z.surv.s <- NULL
    if(nbevt > 0)
      z.surv.s <- surv.s

    
    X0 <- cbind(Xtime, Xfixed.Y, Xfixed.s, Xsurvcause, Xsurv.s)
    
    nom.unique <- unique(colnames(X0))
    X0 <- X0[,nom.unique,drop=FALSE]
    ##X0 fini
    
    
    ##ordonner les mesures par individu
    matYX <- cbind(IND,Y0,indiceY0,outcome,X0)
    matYXord <- matYX[order(IND),]
    Y0 <- as.numeric(matYXord[,2])
    X0 <- apply(matYXord[,-c(1:4),drop=FALSE],2,as.numeric)
    IND <- matYXord[,1]
    outcome <- matYXord[,4]
    indiceY0 <- as.numeric(matYXord[,3])
    
    
    
    if(nbevt>0){
      
      dataSurv <- dataSurv[which(dataSurv[,1] %in% IND),]
      dataSurv <- dataSurv[order(dataSurv[,1]),]
      nmes <- as.vector(table(dataSurv[,1]))
      data.surv <- apply(dataSurv[cumsum(nmes),-1],2,as.numeric) # 1 line per subject
      tsurv0 <- data.surv[,1]
      tsurv <- data.surv[,2]
      devt <- data.surv[,3]
      
      
      ## test de hazard
      arghaz <- hazard
      hazard <- rep(hazard,length.out=nbevt)
      if(any(hazard %in% c("splines","Splines")))
      {
        hazard[which(hazard %in% c("splines","Splines"))] <- "5-quant-splines"
      }
      if(any(hazard %in% c("piecewise","Piecewise")))
      {
        hazard[which(hazard %in% c("piecewise","Piecewise"))] <- "5-quant-piecewise"
      }
      
      haz13 <- strsplit(hazard[which(!(hazard=="Weibull"))],"-")
      if(any(sapply(haz13,length)!=3)) stop("Invalid argument hazard")
      
      nz <- rep(2,nbevt)
      locnodes <- NULL
      typrisq <- rep(2,nbevt)
      nprisq <- rep(2,nbevt)
      
      nznodes <- 0 #longueur de hazardnodes
      ii <- 0
      if(any(hazard!="Weibull"))
      {
        
        for (i in 1:nbevt)
        {
          if(hazard[i]=="Weibull") next;
          
          ii <- ii+1
          
          nz[i] <- as.numeric(haz13[[ii]][1])
          if(nz[i]<3) stop("At least 3 nodes are required")
          typrisq[i] <- ifelse(haz13[[ii]][3] %in% c("splines","Splines"),3,1)
          nprisq[i] <- ifelse(haz13[[ii]][3] %in% c("splines","Splines"),nz[i]+2,nz[i]-1)
          locnodes <- c(locnodes, haz13[[ii]][2])
          if(!(haz13[[ii]][3] %in% c("splines","Splines","piecewise","Piecewise"))) stop("Invalid argument hazard")
          
          if((haz13[[ii]][2]=="manual"))
          {
            nznodes <- nznodes + nz[i]-2
          }
          
          if(!all(locnodes %in% c("equi","quant","manual"))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")
        }
        
        if(!is.null(hazardnodes))
        {
          if(!any(locnodes == "manual"))  stop("hazardnodes should be NULL if the nodes are not chosen manually")
          
          if(length(hazardnodes) != nznodes) stop(paste("Vector hazardnodes should be of length",nznodes))
        }
      }
      else
      {
        if(!is.null(hazardnodes)) stop("hazardnodes should be NULL if Weibull baseline risk functions are chosen")
      }
      
      
      if(nbevt>1 & length(arghaz)==1 & nznodes>0)
      {
        hazardnodes <- rep(hazardnodes,length.out=nznodes*nbevt)
      }
      
      nrisqtot <- sum(nprisq)
      
      zi <- matrix(0,nrow=max(nz),ncol=nbevt)
      nb <- 0
      
      minT1 <- 0
      maxT1 <- max(tsurv)
      tsurvevt <- tsurv
      
      if(idtrunc==1)
      {
        minT1 <- min(tsurv,tsurv0)
        maxT1 <- max(tsurv,tsurv0)
      }
      
      ## arrondir
      minT2 <- round(minT1,3)
      if(minT1<minT2) minT2 <- minT2-0.001
      minT <- minT2
      
      maxT2 <- round(maxT1,3)
      if(maxT1>maxT2) maxT2 <- maxT2+0.001
      maxT <- maxT2
      
      if(length(hazardrange)){
        if(hazardrange[1]>minT) stop(paste("hazardrange[1] should be <=",minT))
        if(hazardrange[2]<maxT) stop(paste("hazardrange[2] should be >=",maxT))
        minT <- hazardrange[1]
        maxT <- hazardrange[2]
      }
      
      startWeib <- rep(0,nbevt)
      startWeib[which(typrisq==2)] <- rep(startWeibull, length.out=length(which(typrisq==2)))
      ii <- 0
      for(i in 1:nbevt)
      {
        if(typrisq[i]==2)
        {
          if(minT < startWeib[i]) stop("Some entry or event times are bellow startWeibull")
          zi[1:2,i] <- c(startWeib[i],maxT)
        }
        else
        {
          ii <- ii+1
          
          if(locnodes[ii]=="manual")
          {
            zi[1:nz[i],i] <- c(minT,hazardnodes[nb+1:(nz[i]-2)],maxT)
            nb <- nb + nz[i]-2
          }
          if(locnodes[ii]=="equi")
          {
            zi[1:nz[i],i] <- seq(minT,maxT,length.out=nz[i])
          }
          if(locnodes[ii]=="quant")
          {
            pi <- c(1:(nz[i]-2))/(nz[i]-1)
            qi <- quantile(tsurvevt,prob=pi)
            zi[1,i] <- minT
            zi[2:(nz[i]-1),i] <- qi
            zi[nz[i],i] <- maxT
          }
        }
      }
      
      
      ## structure of association
      if(length(sharedtype)!=(nbevt*ny)) stop("In sharedtype, one type should be specified per outcome and per event.")
      if(any(!(sharedtype %in% c("RE","CV")))) stop("The value of argument sharedtype must be 'RE' (for shared random effects) or 'CV' (for outcome current value)")
      if(any(sharedtype == "CV") & !(nGK %in% c(7,15))) stop("The value of argument nGK should be rather 7 or 15.")
      
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
    
    
    ##remplir idg, idea, etc
    z.X0 <- strsplit(nom.unique[-1],split=":",fixed=TRUE)
    z.X0 <- lapply(z.X0,sort)
    
    idg <- NULL
    idea <- NULL
    for(k in 1:ny){ # idg : 1 si EF, 2 si EF en interaction ac tps
      idg <- rbind(idg,
                   c(0,(z.X0 %in% z.fixed.Y[[k]] + 0)))
      idea <- c(idea,z.random.Y[k])
    }
    idg.interac <- fixed.Y.interac
    
    idg.s <- c(0,(z.X0 %in% z.fixed.s) + 0)
    
    idsurv <- c(0,(z.X0 %in% z.survcause) + 0)
    idsurv[2] <- 0 # no intercept or confused with baseline risk
    
    idsurv.s <- 0
    if(nbevt > 0)
      idsurv.s <- z.surv.s
      
    # id.fp = fct du tps
    FP.power <- c(-2,-1,-0.5,0,0.5,1,2,3)
    if(!is.null(unlist(FP.power.fx)) & !any(unlist(FP.power.fx) %in% FP.power)) stop(paste("FP powers should be among :",paste(FP.power,collapse=", ")))
    if(!is.null(unlist(FP.power.rd)) & !any(unlist(FP.power.rd) %in% FP.power)) stop(paste("FP powers should be among :",paste(FP.power,collapse=", ")))
    idg.fp <- NULL
    idea.fp <- NULL
    for(k in 1:ny){ # 2 lignes / outcome : EF et EA; 1 si x^p, 2 si x^p + x^p*log(x)
      idg.fp <- rbind(idg.fp,
                      FP.power %in% FP.power.fx[[k]] + 0)
      idea.fp <- rbind(idea.fp,
                       FP.power %in% FP.power.rd[[k]] + 0)
    }
    if(any(idg.fp>2)) stop("In FP() functions, one power cannot appear more than 2 times.")
    if(any(idea.fp>2)) stop("In FP() functions, one power cannot appear more than 2 times.")
    
    idsurv.s.fp <- NULL
    if(!is.null(FP.power.s) & !any(FP.power.s %in% FP.power)) stop(paste("FP powers should be among :",paste(FP.power,collapse=",")))
    idsurv.s.fp <- FP.power %in% FP.power.s + 0
    if(any(idsurv.s.fp>2)) stop("In FP() functions, one power cannot appear more than 2 times.")
    
    
    
    ## nb prms
    
    # coef pred lin pour survie
    nvarbyevt <- sum(idsurv==1)
    nvarxevt <- nbevt*nvarbyevt #sum(idsurv==1) + nbevt*sum(idsurv==2)
    
    # partie long
    nef.fp <- rowSums(idg.fp)
    nef.var <- rowSums(idg)
    nef.interac <- idg.interac * nef.fp * nef.var
    nef <- nef.fp + nef.var + nef.interac
    nef <- nef + ifelse(idmodel==0,1,0) # logistic +1 pr nu
    
    nea.fp <- rowSums(idea.fp)
    nea.intercept <- idea
    nea <- nea.fp + nea.intercept
    nea <- nea + ifelse(idmodel==0,2,0) # logistic +2 pr nu et rate
    
    nvc <- ifelse(idmodel==1,nea*(nea+1)/2,nea) # logistic : EA indep - linear : EA corr
    nvc <- nvc + ifelse(idmodel==1 & idlink!=0,-1,0) # linear + link fct : 1st RE variance fixed to 1 (identifiability constraint)
    
    ntr <- rep(0,ny)
    ntr[which(idlink==1)] <- 2 # linear
    ntr[which(idlink==2)] <- nbzitr[which(idlink==2)]+2 # spl
    ntrtot <- sum(ntr)
    
    nvc.err <- ny # erreur
    
    # shift
    nef.s <- sum(idg.s==1) # time shift
    nvc.s <- 1
    
    # nb EAs
    neatot <- sum(nea) + nvc.s # nea par mq, +1 pr s
    
    ## nb parametres d'association
    nasso <- 0
    sharedtype <- ifelse(sharedtype=="RE",0,1)
    if(nbevt>0){
      
      nasso <- rep(0,nbevt)
      
      # markers' characteristics 
      sharedtype <- matrix(sharedtype,nrow=nbevt,ncol=ny,byrow=TRUE)
      
      for(ke in 1:nbevt){
        for(k in 1:ny){
          if(sharedtype[ke,k]==0) # REs
            nasso[ke] <- nasso[ke] + nea[k]
          if(sharedtype[ke,k]==1) # CV
            nasso[ke] <- nasso[ke] + 1
        }
      }
      
      # shift
      nasso.s <- surv.s
      nasso.s.fp <- sum(idsurv.s.fp)
      
      nasso <- nasso + nasso.s + nasso.s.fp
      
    }
    
    ##nombre total de parametres
    NPM <- nrisqtot + nvarxevt + sum(nasso) +
      nef.s + nvc.s +
      sum(nef) + sum(nvc) + ntrtot + nvc.err
    
    V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres
    
    
    ## parametres MC
    methInteg <- switch(methInteg,"MCO"=1,"MCA"=2,"QMC"=3)
    seqMC <- 0
    dimMC <- 0
    if(methInteg==3)
    {
      dimMC <- neatot
      if(dimMC>0) seqMC <- randtoolbox::sobol(n=nMC,dim=dimMC,normal=TRUE,scrambling=1)
    }
    #print("seqMC");cat(seqMC)
    
    
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
      
      if(nbevt>0){  #si Weibull et prms au carre pr positivite -> valeurs par defaut = 1
        if(any(hazard!="Weibull")==FALSE & isFALSE(logscale)==TRUE){
          sumrisq = 0
          for(i in 1:nbevt){
            if(typrisq[i]==2){
              b[sumrisq+1:nprisq[i]] <- 1
              if(idtrunc==1 & any(Tentry==0)) #si entree retardee et au moins un temps vaut 0
                b[sumrisq+nprisq[i]] <- 1.25 #sinon pblm lors de recherche prms, risq instant tend vers infini qd w2 < 1 et t=0
            }
            sumrisq = sumrisq + nprisq[i] + nvarbyevt + nasso[i]
          }
        }
        if(any(hazard %in% c("splines","Splines")) & idtrunc==1 & any(Tentry==0)){
          sumrisq = 0
          for(i in 1:nbevt){
            if(typrisq[i]==3){
              b[sumrisq+1:nprisq[i]] <- 10^-7 #sinon pgrm crash
            }
            sumrisq = sumrisq + nprisq[i] + nvarbyevt + nasso[i]
          }
        }
      }
      
      #variances to 1
      b[nrisqtot + nvarxevt + sum(nasso) + nef.s + nvc.s] <- 1 # shift
      for(k in 1:ny){
        prev_k <- ifelse(k==1,0,1)*(sum(nef[1:(k-1)])+sum(nvc[1:(k-1)])+sum(ntr[1:(k-1)])+(k-1)) # previous markers' prms
        if(idmodel[k]==0){ # logistic
          b[nrisqtot + nvarxevt + sum(nasso) + nef.s + nvc.s + prev_k + nef[k] + 1:nea[k]] <- 1 # var rate, nu et u
        }
        if(idmodel[k]==1 & nea[k]>1){ # linear
          init.nvc <- diag(nea[k])
          init.nvc <- init.nvc[upper.tri(init.nvc, diag=TRUE)]
          if(idlink[k]!=0) init.nvc <- init.nvc[-1]
          b[nrisqtot + nvarxevt + sum(nasso) + nef.s + nvc.s + prev_k + nef[k] + 1:nvc[k]] <- init.nvc
        }
        b[nrisqtot + nvarxevt + sum(nasso) + nef.s + nvc.s + prev_k + nef[k] + nvc[k] + ntr[k] + 1] <- 1
      }
      
      #fct de lien
      if(any(idlink!=0)){
        for(k in 1:ny){
          prev_k <- ifelse(k==1,0,1)*(sum(nef[1:(k-1)])+sum(nvc[1:(k-1)])+sum(ntr[1:(k-1)])+(k-1)) # previous markers' prms
          if(idlink[k]==1){ # linear
            b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+nvc[k]+1] <- mean(get(dataY[k])[,nomsY[k]])
            b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+nvc[k]+2] <- 1
          }
          if(idlink[k]==2){ # spl
            b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+nvc[k]+1] <- -2
            b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+nvc[k]+2:ntr[k]] <- 0.1
          } 
        }
      }
    }
    
    
    ##------------------------------------------
    ##------nom au vecteur best
    ##--------------------------------------------
    
    nom.X0 <- colnames(X0)[-1]
    nom.X0[nom.X0=="(Intercept)"] <- "intercept"
    
    if(nbevt>0)
    {
      
      ##prm fct de risque
      if(isTRUE(logscale))
      {
        sumrisq = 0
        for(i in 1:nbevt)
        {
          nom1 <- rep(paste("event",i,sep=""),nprisq[i])
          if(typrisq[i]==2)
          {
            names(b)[sumrisq+1:nprisq[i]] <- paste(nom1[1:2]," log(Weibull",1:2,")",sep="")
          }
          if(typrisq[i]==1)
          {
            names(b)[sumrisq+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(piecewise",1:(nz[i]-1),")",sep="")
          }
          if(typrisq[i]==3)
          {
            names(b)[sumrisq+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(splines",1:(nz[i]+2),")",sep="")
          }
          sumrisq = sumrisq + nprisq[i] + nvarbyevt + nasso[i]
        }
      }
      else
      {
        sumrisq = 0
        for(i in 1:nbevt)
        {
          nom1 <- rep(paste("event",i,sep=""),nprisq[i])
          if(typrisq[i]==2)
          {
            names(b)[sumrisq+1:nprisq[i]] <- paste(nom1[1:2]," +/-sqrt(Weibull",1:2,")",sep="")
          }
          if(typrisq[i]==1)
          {
            names(b)[sumrisq+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(piecewise",1:(nz[i]-1),")",sep="")
          }
          if(typrisq[i]==3)
          {
            names(b)[sumrisq+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(splines",1:(nz[i]+2),")",sep="")
          }
          sumrisq = sumrisq + nprisq[i] + nvarbyevt + nasso[i]
        }
      }
      
      
      ##prm covariables survival
      nom1 <- NULL
      for(j in 1:nv)
      {
        if(idsurv[j]==1) # pred lin
        {
          nom1 <- c(nom1,c("NA",nom.X0)[j])
        }
      }
      if(nvarxevt>0){ 
        sumrisq = 0
        for(i in 1:nbevt)
        {
          names(b)[sumrisq+nprisq[i]+1:nvarbyevt] <- paste("event",i," ",nom1,sep="")
          sumrisq = sumrisq + nprisq[i] + nvarbyevt + nasso[i]
        }
      }

      ##prm association
      sumrisq = 0
      for(i in 1:nbevt){
        
        # shift
        if(nasso.s)
          names(b)[sumrisq+nprisq[i]+nvarbyevt+nasso.s] <- paste("event",i," asso.s intercept",sep="")
        if(nasso.s.fp){
          # time fct
          fp <- NULL
          for(l in 1:length(FP.power)){
            if(idsurv.s.fp[l]!=0){
              fp <- c(fp,as.character(FP.power[l]))
              if(idsurv.s.fp[l]==2)
                fp <- c(fp,paste(FP.power[l],".log",sep=""))
            }
          }
          names(b)[sumrisq+nprisq[i]+nvarbyevt+nasso.s+1:nasso.s.fp] <- paste("event",i," asso.s FP[",fp,"]",sep="")
        }
          
        
        # markers' characteristics
        
        sumasso <- 0
        
        for(k in 1:ny){
          
          if(sharedtype[i,k] == 0){ # RE
            if(idmodel[k] == 0) # logistic
              names(b)[sumrisq+nprisq[i]+nvarbyevt+nasso.s+nasso.s.fp+sumasso+1:nea[k]] <- paste("event",i," asso.Y",k," re.",c("rate","nu","int"),sep="")
            if(idmodel[k] == 1) # linear
              names(b)[sumrisq+nprisq[i]+nvarbyevt+nasso.s+nasso.s.fp+sumasso+1:nea[k]] <- paste("event",i," asso.Y",k," re",1:nea[k],sep="") 
            sumasso <- sumasso + nea[k]
          }
          
          if(sharedtype[i,k] == 1){ # CV
            names(b)[sumrisq+nprisq[i]+nvarbyevt+nasso.s+nasso.s.fp+sumasso+1] <- paste("event",i," asso.Y",k," cv",sep="")
            sumasso <- sumasso + 1
          }
          
        }
        sumrisq = sumrisq + nprisq[i] + nvarbyevt + nasso[i]
      }
      
    }
    
    # shift
    names(b)[nrisqtot+nvarxevt+sum(nasso)+1:nef.s] <- paste("s. ",nom.X0[idg.s[-1]!=0],sep="")
    names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+1] <- "s. var"
    
    
    # partie long
    for(k in 1:ny){
      
      prev_k <- ifelse(k==1,0,1)*(sum(nef[1:(k-1)])+sum(nvc[1:(k-1)])+sum(ntr[1:(k-1)])+(k-1)) # previous markers' prms
      
      
      if(idmodel[k]==0){ # logistic
        names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+1:(nef[k]-1)] <- paste("Y",k,". rate. ",nom.X0[idg[k,-1]!=0],sep="")
        names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]] <- paste("Y",k,". v. mu",sep="")
        names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nvc[k]] <- paste("Y",k,". ",c("rate","v","u"),". var",sep="")
        names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+nvc[k]+1] <- paste("Y",k,". std.err",sep="")
      }
      
      if(idmodel[k]==1){ # linear
        
        ## fixed effects
        # time fct
        fp <- NULL
        for(l in 1:length(FP.power)){
          if(idg.fp[k,l]!=0){
            fp <- c(fp,as.character(FP.power[l]))
            if(idg.fp[k,l]==2)
              fp <- c(fp,paste(FP.power[l],".log",sep=""))
          }
        }
        names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+1:nef.fp[k]] <- paste("Y",k,". FP[",fp,"]",sep="")
        # covariates
        if(nef.var[k])
          names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef.fp[k]+1:nef.var[k]] <- paste("Y",k,". ",nom.X0[idg[k,-1]!=0],sep="")
        # interaction
        if(nef.interac[k])
          names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef.fp[k]+nef.var[k]+1:nef.interac[k]] <- paste("Y",k,". FP[",paste(rep(fp, each = length(nom.X0[idg[k,-1]!=0])),"]:", nom.X0[idg[k,-1]!=0], sep = ""),sep="")
        
        names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nvc[k]] <- paste("Y",k,". varcov ",1:nvc[k],sep="")
        names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+nvc[k]+ntr[k]+1] <- paste("Y",k,". std.err",sep="")
        
      } 
        
      if(idlink[k]==1) # linear
        names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+nvc[k]+1:ntr[k]] <- paste("Y",k,". Linear",1:ntr[k],sep="")
      if(idlink[k]==2) # spl
        names(b)[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+nvc[k]+1:ntr[k]] <- paste("Y",k,". I-splines",1:ntr[k],sep="")
      
    }
    
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
    
    
    ## RE : varcov -> cholesky (= std si diago / pas varcov)
    mvc <- matrix(0,neatot,neatot)
      #shift
    mvc[1,1] <- b[nrisqtot + nvarxevt + sum(nasso) + nef.s + nvc.s]
      #longit.
    for(k in 1:ny){
      prev_k <- ifelse(k==1,0,1)*(sum(nef[1:(k-1)])+sum(nvc[1:(k-1)])+sum(ntr[1:(k-1)])+(k-1)) # previous markers' prms
      if(idmodel[k]==0){ # logistic
        diag(mvc)[1+ifelse(k==1,0,1)*sum(nea[1:(k-1)])+1:nea[k]] <- b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nea[k]]
      }
      if(idmodel[k]==1){ # linear
        mvc.k <- matrix(0,nea[k],nea[k])
        b_mvc.k <- NULL
        if((idlink[k]==0) | (idlink[k]!=0 & nea[k]>1)) b_mvc.k <- b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nvc[k]]
        if(idlink[k]!=0) b_mvc.k <- c(1,b_mvc.k) # identifiability constraint
        mvc.k[upper.tri(mvc.k,diag=T)] <- b_mvc.k
        mvc.k <- t(mvc.k)
        mvc.k[upper.tri(mvc.k, diag=TRUE)] <- b_mvc.k
        mvc[1+ifelse(k==1,0,1)*sum(nea[1:(k-1)])+1:nea[k],1+ifelse(k==1,0,1)*sum(nea[1:(k-1)])+1:nea[k]] <- mvc.k
      }
    }
    ind_var0 <- (1:neatot)[diag(mvc)==0] # si variances a 0 -> les passer en 1 pr transfo cholesky [uniquement valable pr shift ou logistic car diago]
    if(length(ind_var0)) diag(mvc)[ind_var0] <- 1 # variances 0 -> 1 pour la transfo cholesky
    ch <- chol(mvc)
    if(length(ind_var0)) diag(ch)[ind_var0] <- 0 # variances 1 -> 0 apres la transfo cholesky
    # re-injecter
    b[nrisqtot + nvarxevt + sum(nasso) + nef.s + nvc.s] <- ch[1,1]
    for(k in 1:ny){
      prev_k <- ifelse(k==1,0,1)*(sum(nef[1:(k-1)])+sum(nvc[1:(k-1)])+sum(ntr[1:(k-1)])+(k-1)) # previous markers' prms
      if(idmodel[k]==0) # logistic
        b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nea[k]] <- diag(ch)[1+ifelse(k==1,0,1)*sum(nea[1:(k-1)])+1:nea[k]]
      if(idmodel[k]==1){ # linear
        ch.k <- ch[1+ifelse(k==1,0,1)*sum(nea[1:(k-1)])+1:nea[k],1+ifelse(k==1,0,1)*sum(nea[1:(k-1)])+1:nea[k]]
        ch.k_vect <- ch.k[upper.tri(ch.k,diag=T)]
        if(idlink[k]!=0) ch.k_vect <- ch.k_vect[-1] # identifiability constraint
        b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nvc[k]] <- ch.k_vect
      }
    }
    
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
    
    
    ######################################################################################
    ######################################################################################
    ######################################################################################
    
    ###############
    ###   MLA   ###
    ###############
    
    # print("b");print(b)
    # print("Y0");print(Y0)
    # print("X0");print(X0)
    # print("tsurv0");print(tsurv0)
    # print("tsurv");print(tsurv)
    # print("devt");print(devt)
    # print("idea.fp");print(idea.fp)
    # print("idg.fp");print(idg.fp)
    # print("idea");print(idea)
    # print("idg");print(idg)
    # print("idg.interac");print(idg.interac)
    # print("idg.s");print(idg.s)
    # print("idsurv");print(idsurv)
    # print("idsurv.s");print(idsurv.s)
    # print("idsurv.s.fp");print(idsurv.s.fp)
    # print("sharedtype");print(sharedtype)
    # print("nGK");print(nGK)
    # print("typrisq");print(typrisq)
    # print("nz");print(nz)
    # print("zi");print(zi)
    # print("nbevt");print(nbevt)
    # print("idtrunc");print(idtrunc)
    # print("logspecif");print(logspecif)
    # print("ny");print(ny)
    # print("ns");print(ns)
    # print("nv");print(nv)
    # print("nobs");print(nobs)
    # print("nmes");print(nmes)
    # print("NPM");print(NPM)
    # print("nfix");print(nfix)
    # print("bfix");print(bfix)
    # print("idmodel");print(idmodel)
    # print("idlink");print(idlink)
    # print("nbzitr");print(nbzitr)
    # print("zitr");print(zitr)
    # print("uniqueY0");print(uniqueY0)
    # print("indiceY0");print(indiceY0)
    # print("nvalSPLORD");print(nvalSPLORD)
    # print("fix");print(fix)
    # print("methInteg");print(methInteg)
    # print("nMC");print(nMC)
    # print("dimMC");print(dimMC)
    # print("seqMC");print(seqMC)
    
    
    if(maxiter==0)
    {
      vrais <- loglik(b,Y0,X0,tsurv0,tsurv,devt,
                      idea.fp,idg.fp,idea,idg,idg.interac,idg.s,idsurv,idsurv.s,idsurv.s.fp,sharedtype,nGK,
                      typrisq,nz,zi,nbevt,idtrunc,logspecif,
                      ny,ns,nv,
                      nobs,nmes,NPM,nfix,bfix,
                      idmodel,idlink,nbzitr,zitr,uniqueY0,indiceY0,nvalSPLORD,fix,
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
                 Y0=Y0,X0=X0,Tentr0=tsurv0,Tevt0=tsurv,Devt0=devt,
                 ideafp0=idea.fp,idgfp0=idg.fp,idea0=idea,idg0=idg,idginterac0=idg.interac,idgs0=idg.s,idsurv0=idsurv,idsurvs0=idsurv.s,idsurvsfp0=idsurv.s.fp,sharedtype0=sharedtype,nGK0=nGK,
                 typrisq0=typrisq,nz0=nz,zi0=zi,nbevt0=nbevt,idtrunc0=idtrunc,logspecif0=logspecif,
                 ny0=ny,ns0=ns,nv0=nv,
                 nobs0=nobs,nmes0=nmes,npm0=NPM,nfix0=nfix,bfix0=bfix,
                 idmodel0=idmodel,idlink0=idlink,nbzitr0=nbzitr,zitr0=zitr,uniqueY0=uniqueY0,indiceY0=indiceY0,nvalSPLORD0=nvalSPLORD,fix0=fix,
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
    
    
    ## RE : cholesky or std -> varcov 
    Cholesky <- NULL
    #shift
    Cholesky <- c(Cholesky,out$best[nrisqtot + nvarxevt + sum(nasso) + nef.s + nvc.s])
    out$best[nrisqtot + nvarxevt + sum(nasso) + nef.s + nvc.s] <- out$best[nrisqtot + nvarxevt + sum(nasso) + nef.s + nvc.s]**2
    #longit.
    for(k in 1:ny){
      prev_k <- ifelse(k==1,0,1)*(sum(nef[1:(k-1)])+sum(nvc[1:(k-1)])+sum(ntr[1:(k-1)])+(k-1))
      if(idmodel[k]==0){ # logistic
        Cholesky <- c(Cholesky,out$best[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nea[k]])
        out$best[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nea[k]] <- out$best[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nea[k]]**2
      }
      if(idmodel[k]==1){ # linear
        if(nea[k]==1){
          if(idlink[k]==0){
            Cholesky <- c(Cholesky,out$best[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1])
            out$best[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1] <- out$best[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1]**2
          }
        }
        if(nea[k]>1){
          chol.k_vect <- out$best[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nvc[k]]
          if(idlink[k]!=0) chol.k_vect <- c(1,chol.k_vect)
          Cholesky <- c(Cholesky,chol.k_vect)
          chol.k <- matrix(0,nea[k],nea[k])
          chol.k[upper.tri(chol.k,diag=T)] <- chol.k_vect
          vc.k <- t(chol.k) %*% chol.k
          vc.k_vect <- vc.k[upper.tri(vc.k,diag=T)]
          if(idlink[k]!=0) vc.k_vect <- vc.k_vect[-1]
          out$best[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+prev_k+nef[k]+1:nvc[k]] <- vc.k_vect
        }
      }
    }
    
    ##predictions
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
    
    
    ##estimlink
    ysim <- matrix(out$marker,nsim,ny)
    transfo <- matrix(out$transfY,nsim,ny)
    estimlink <- as.vector(rbind(ysim,transfo))
    estimlink <- matrix(estimlink,nsim,2*ny)
    colnames(estimlink) <- paste(c("","transf"),rep(nomsY, each=2),sep="")
    
    
    # N <- rep(NA,14) # /!\ a revoir au moment de faire le summary
    # N[1] <- ny
    # N[2] <- nobs
    # N[3] <- nbevt
    # N[4] <- nrisqtot
    # N[5] <- nvarxevt
    # N[6] <- sum(nasso)
    # N[7] <- nef.s
    # N[8] <- nvc.s
    # N[9] <- sum(nef)
    # N[10] <- sum(nea)
    # N[11] <- sum(nvc)
    # N[12] <- sum(ntr)
    # N[13] <- nvc.err
    
    N <- list()
    N[[1]] <- ny
    N[[2]] <- colSums(nmes) #nobs
    N[[3]] <- nbevt
    N[[4]] <- nrisqtot
    N[[5]] <- nvarxevt
    N[[6]] <- sum(nasso)
    N[[7]] <- nef.s
    N[[8]] <- nvc.s
    N[[9]] <- nef
    N[[10]] <- nea
    N[[11]] <- nvc
    N[[12]] <- ntr
    N[[13]] <- nvc.err
    
    
    nevent <- rep(0,nbevt)
    for(ke in 1:nbevt)
    {
      nevent[ke] <- length(which(devt==ke))
    }
    
    Nprm <- c(nprisq)
    
    
    nom.X0[nom.X0=="intercept"] <- "Intercept"
    
    
    
    ## noms des variables
    Names <- list(Xnames=nom.X0,Ynames=nomsY,
                  ID=nom.subject,Tnames=noms.surv,
                  Xvar=setdiff(ttesLesVar,noms.surv))
    
    form <- list(link=link, model=model,
                 fixed.s=fixed.s,
                 form.cause=form.cause)
    
    cost <- proc.time()-ptm
    
    
    res <-list(ns=ns,idg=idg,idea=idea,idgfp=idg.fp,ideafp=idea.fp,idginterac=idg.interac,idsurv=idsurv,idsurvs=idsurv.s,idsurvsfp=idsurv.s.fp,
               loglik=out$loglik,best=out$best,V=V,gconv=out$gconv,conv=out$conv,
               call=cl,niter=out$niter,N=N,nevent=nevent,Nprm=Nprm,
               pred=pred,Names=Names,form=form,cholesky=Cholesky,
               logspecif=logspecif,predSurv=predSurv,typrisq=typrisq,hazardnodes=zi,nz=nz, sharedtype=sharedtype,
               estimlink=estimlink,modeltype=idmodel,linktype=idlink,linknodes=zitr,nbnodes=nbnodes,nbzitr=nbzitr,
               na.action=nayk,AIC=2*(length(out$best)-length(posfix)-out$loglik),BIC=(length(out$best)-length(posfix))*log(ns)-2*out$loglik,data=datareturn,
               #wRandom=wRandom,b0Random=b0Random,
               posfix=posfix,CPUtime=cost[3])
    
    names(res$best) <- namesb
    class(res) <-c("jointLTSM")
    
    return(res)
    
    
}


loglik <- function(b0,Y0,X0,Tentr0,Tevt0,Devt0,
                   ideafp0,idgfp0,idea0,idg0,idginterac0,idgs0,idsurv0,idsurvs0,idsurvsfp0,sharedtype0,nGK0,
                   typrisq0,nz0,zi0,nbevt0,idtrunc0,logspecif0,
                   ny0,ns0,nv0,
                   nobs0,nmes0,npm0,nfix0,bfix0,
                   idmodel0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0,nvalSPLORD0,fix0,
                   methInteg0,nMC0,dimMC0,seqMC0)
{
    res <- 0
    .Fortran(C_loglik,as.double(Y0),as.double(X0),as.double(Tentr0),as.double(Tevt0),as.integer(Devt0),
             as.integer(ideafp0),as.integer(idgfp0),as.integer(idea0),as.integer(idg0),as.integer(idginterac0),as.integer(idgs0),as.integer(idsurv0),as.integer(idsurvs0),as.integer(idsurvsfp0),as.integer(sharedtype0),as.integer(nGK0),
             as.integer(typrisq0),as.integer(nz0),as.double(zi0),as.integer(nbevt0),as.integer(idtrunc0),as.integer(logspecif0),
             as.integer(ny0),as.integer(ns0),as.integer(nv0),
             as.integer(nobs0),as.integer(nmes0),as.integer(npm0),as.double(b0),as.integer(nfix0),as.double(bfix0),
             as.integer(idmodel0),as.integer(idlink0),as.integer(nbzitr0),as.double(zitr0),as.double(uniqueY0),as.integer(indiceY0),as.integer(nvalSPLORD0),as.integer(fix0),
             as.integer(methInteg0),as.integer(nMC0),as.integer(dimMC0),as.double(seqMC0),loglik_res=as.double(res))$loglik_res
}
