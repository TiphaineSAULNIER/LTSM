#' Predictions of the random-effects
#' 
#' The function computes the predicted values of the random effects given observed data
#' provided in input. 
#' 
#' @param mod an object inheriting from class \code{jointLTSM}.
#' @param newdata data frame containing the data from which predictions are to be computed. 
#' The data frame should include at least all the covariates and 
#' the marker(s) values and the grouping structure. Names should match exactly the names 
#' of the variables in the model.
#' @param subject character specifying the name of the grouping structure.
#' If NULL (the default), the same as in the model will be used.
#' @param re.init optional vector of random effects, default is a random pick according to 
#' random-effects distribution  
#' @param re.posfix optional vector giving the indices in vector re.init of the
#' random effects that should not be estimated. Default to NULL, all random effects are
#' estimated.
#' 
#' @return a list containing:
#' \item{ns}{number of grouping units in the dataset}
#' \item{best}{vector of parameter estimates}
#' \item{predRE}{table containing individual predictions of the random-effects : 
#' a line per subject.}
#' \item{predRE_varcov}{list of vectors containing the upper triangle matrix of variance-covariance
#' estimates of the individually predicted random effects}
#' 
#' @author Tiphaine Saulnier, Cecile Proust-Lima 
#' 
#' @export 


predictRE.LTSM <- function(mod, newdata, subject=NULL, re.init=NULL, re.posfix=NULL){
  
  if(missing(mod)){ stop("The argument mod should be specified")}
  if(class(mod)!="jointLTSM") stop("The argument mod must be an element inheriting from class 'jointLTSM'")
  if(missing(newdata)){ stop("The argument newdata should be specified")}
  if(nrow(newdata)==0) stop("newdata should not be empty")
  
  if(is.null(subject)) subject <- mod$call$subject
  if(!(subject %in% colnames(newdata))) stop("Unable to find variable 'subject' in 'newdata'")
  if(!is.numeric(newdata[,subject])) stop("The argument subject must be numeric in newdata")
  nom.subject <- as.character(subject)
  
  fixed <- eval(mod$call$fixed)
  random <- eval(mod$call$random)
  nlp <- length(fixed) # nb of latent processes
  
  
  link <- mod$form$link
  model <- mod$form$model
  if(!any(model %in% c("logistic","linear"))) stop("The argument model must only contain 'logistic' or 'linear'")
  timevar <- mod$call$timevar
  if(!(timevar %in% colnames(newdata))) stop("Unable to find variable 'timevar' in 'newdata'")
  fixed.s <- mod$form$fixed.s
  
  data <- newdata
  
  B <- mod$best
  
  # verification convergence
  if(mod$conv!=1) stop("The model should have converged to predict.")
  
  
  # arguments jointLTSM
  # si non-renseigne : valeurs par defaut
  
  range <- rep(list(NULL),nlp); if(!is.null(mod$call$range)) range <- eval(mod$call$range)
  intnodes <- NULL; if(!is.null(mod$call$intnodes)) intnodes <- eval(mod$call$intnodes)
  
  na.action <- ifelse(!is.null(mod$call$na.action),mod$call$na.action,1)
  
  survival <- NULL; if(!is.null(mod$call$survival)){ survival <- eval(mod$call$survival); class(survival) <- "formula"}
  #hazard <- ifelse(!is.null(mod$call$hazard),eval(mod$call$hazard),"Weibull")
  hazard <- "Weibull"; if(!is.null(mod$call$hazard)) hazard <- eval(mod$call$hazard)
  hazardrange <- NULL; if(!is.null(mod$call$hazardrange)) hazardrange <- eval(mod$call$hazardrange)
  hazardnodes <- NULL; if(!is.null(mod$call$hazardnodes)) hazardnodes <- eval(mod$call$hazardnodes)
  logscale <- ifelse(!is.null(mod$call$logscale),eval(mod$call$logscale),FALSE)
  startWeibull <- ifelse(!is.null(mod$call$startWeibull),eval(mod$call$startWeibull),0)
  surv.shift <- NULL; if(!is.null(mod$call$surv.shift)) surv.shift <- eval(mod$call$surv.shift)
  sharedtype <- "RE"; if(!is.null(mod$call$sharedtype)) sharedtype <- eval(mod$call$sharedtype)
  nGK <- ifelse(!is.null(mod$call$nGK),mod$call$nGK,15)
  
  methInteg <- ifelse(!is.null(mod$call$methInteg),eval(mod$call$methInteg),"QMC")
  nMC <- ifelse(!is.null(mod$call$nMC),eval(mod$call$nMC),1000)
  maxiter <- ifelse(!is.null(mod$call$maxiter),eval(mod$call$maxiter),100)
  if(maxiter==0) stop("Argument maxiter should be greater than 0.")
  convB <- ifelse(!is.null(mod$call$convB),eval(mod$call$convB),0.0001)
  convL <- ifelse(!is.null(mod$call$convL),eval(mod$call$convL),0.0001)
  convG <- ifelse(!is.null(mod$call$convG),eval(mod$call$convG),0.0001)
  partialH <- NULL; if(!is.null(mod$call$partialH)) partialH <- eval(mod$call$partialH)
  nsim <- ifelse(!is.null(mod$call$nsim),eval(mod$call$nsim),100)
  
  verbose <- ifelse(!is.null(mod$call$verbose),eval(mod$call$verbose),FALSE)
  #nproc <- ifelse(!is.null(mod$call$nproc),eval(mod$call$nproc),1)
  
  ##############################################################################
  
  #############################################
  ###   preparation des arguments Fortran   ###
  #############################################
  # strictement identique au code jointLTSM
  
  
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
    surv <- survival[[2]]
    
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
    # nbevt <- length(attr(surv,"states"))
    # if(nbevt<1) stop("No observed event in the data")
    nbevt <- mod$N[[3]]
    
    
    ## pour la formule pour survival
    form.surv <- survival[3]
    
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
  
  
  ##remplir zitr
  zitr <- mod$linknodes
  
  
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
  
  X0 <- cbind(Xtime, Xfixed.Y, Xrandom.Y, Xfixed.s, Xsurvcause, Xsurv.s)
  
  colnames(X0)[1] <- timevar
  nom.unique <- c(timevar,unique(colnames(X0)[-1]))
  X0 <- X0[,nom.unique,drop=FALSE]
  X0 <- as.matrix(X0)
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
    if(length(nmes)==1)
      data.surv <- array(data.surv,dim=c(1,3))
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
  nef.s <- sum(idg.s==1) #time shift
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
  
  ###valeurs des parametres 
  if (length(B)==NPM) b <- B
  else stop(paste("Vector B should be of length",NPM))
  namesb <- names(b)
  
  
  ## RE : varcov
  # pas de transformation ici
  
  ##############################################################################
  
  ## prm fixes
  fix <- rep(0,neatot)
  if(length(re.posfix))
  {
    if(any(!(re.posfix %in% 1:neatot))) stop("Indexes in posfix are not correct")
    if(is.null(re.init)) stop("Argument re.init should be specified if re.posfix is")
    
    fix[re.posfix] <- 1
  }
  if(length(re.posfix)==neatot) stop("No parameter to estimate")
  
  ##############################################################################
  
  # valeurs initiales en fct des ecart-types des prms estimes
  #re <- rep(0.1,neatot)
  if(!is.null(re.init)){
    if(length(re.init)!=neatot) stop(paste("Argument re.init should be of length",neatot))
    re <- re.init
  }
  else{
    # tirage aleatoire selon N(0,V)
    V <- matrix(0,neatot,neatot)
    V[1,1] <- b[nrisqtot+nvarxevt+sum(nasso)+nef.s+1]
    sumPrmK = 0
    sumRE = 1
    for(k in 1:ny){
      
      if(idmodel[k]==0) #logistic
        diag(V[sumRE+1:nea[k],sumRE+1:nea[k]]) <- b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+sumPrmK+nef[k]+1:nea[k]]
      
      if(idmodel[k]==1){ #linear
        U <- matrix(0,nea[k],nea[k])
        U_vect <- b[nrisqtot+nvarxevt+sum(nasso)+nef.s+nvc.s+sumPrmK+nef[k]+1:nvc[k]]
        if(idlink[k]!=0) U_vect <- c(1,U_vect) #linear or spl : identifiability constraint
        U[upper.tri(U,diag=TRUE)] <- U_vect
        U[lower.tri(U,diag=TRUE)] <- U_vect
        V[sumRE+1:nea[k],sumRE+1:nea[k]] <- U
      }
      
      #incrementation
      sumPrmK = sumPrmK + nef[k] + nvc[k] + ntr[k] + 1
      sumRE = sumRE + nea[k]
    }
    
    re <- rmvnorm(n=1,mean=rep(0,neatot),sigma=V); re <- re[1,]
  }  
  
  
  # fixed prms
  nfix <- sum(fix)
  refix <- 0
  if(nfix>0)
  {
    refix <- re[which(fix==1)]
    re <- re[which(fix==0)]
    neatot <- neatot-nfix
  }
  
  ##############################################################################
  
  ##predictions
  predRE <- matrix(NA,nrow=ns,ncol=neatot)
  predRE <- data.frame(unique(IND),rep(NA,ns),rep(NA,ns),predRE)
  predRE.colnames <- "re.s"
  for(k in 1:ny)
    predRE.colnames <- c(predRE.colnames,paste("re.Y",k,".ea",1:nea[k],sep=""))
  colnames(predRE) <- c(nom.subject,"log.dens.value","conv",predRE.colnames[fix==0])
  
  predRE_varcov <- vector("list", ns)
  names(predRE_varcov) <- unique(IND)
  
  ##############################################################################
  
  ###############
  ###   MLA   ###
  ###############
  
  nmescur <- 0
  for(i in 1:ns){
    
    print("#########################")
    print(paste("SUBJECT ",i))
    
    res <- mla(b = re, m = neatot, fn = predREfct,
               clustertype=NULL,.packages=NULL,
               epsa=convB,epsb=convL,epsd=convG,
               digits=8,print.info=verbose,blinding=FALSE,
               multipleTry=25,file="",partialH=partialH,
               nproc=1,maxiter=maxiter,minimize=FALSE,
               best0=b,
               Y0=Y0,X0=X0,Tentr0=tsurv0,Tevt0=tsurv,Devt0=devt,
               ideafp0=idea.fp,idgfp0=idg.fp,idea0=idea,idg0=idg,idginterac0=idg.interac,idgs0=idg.s,idsurv0=idsurv,idsurvs0=idsurv.s,idsurvsfp0=idsurv.s.fp,sharedtype0=sharedtype,nGK0=nGK,
               typrisq0=typrisq,nz0=nz,zi0=zi,nbevt0=nbevt,idtrunc0=idtrunc,logspecif0=logspecif,
               ny0=ny,ns0=ns,nv0=nv,
               nobs0=nobs,nmes0=nmes,npm0=NPM,nea0=neatot,nfix0=nfix,refix0=refix,fix0=fix,
               idmodel0=idmodel,idlink0=idlink,nbzitr0=nbzitr,zitr0=zitr,uniqueY0=uniqueY0,indiceY0=indiceY0,nvalSPLORD0=nvalSPLORD,
               ind0=i,nmescur0=nmescur)
    
    
    predRE[i,2] <- res$fn.value # fct density value
    predRE[i,3] <- res$istop # conv
    predRE[i,-c(1:3)] <- res$b # predRE
    
    predRE_varcov[[i]] <- res$v
    
    
    #incrementation
    nmescur <- nmescur + sum(nmes[i,])
    
    #print(paste("dens value ",res$fn.value))
    #print(paste("conv ",res$istop))
    #print("re ");print(paste(res$b))
    #print(paste("nmescur_i ",nmescur))
    
  }
  
  
  ##############################################################################
  
  res <-list(ns=ns, best=B, predRE=predRE, predRE_varcov=predRE_varcov)
  
  return(res)
  
} 



predREfct <- function(re0,
                      best0,
                      Y0,X0,Tentr0,Tevt0,Devt0,
                      ideafp0,idgfp0,idea0,idg0,idginterac0,idgs0,idsurv0,idsurvs0,idsurvsfp0,sharedtype0,nGK0,
                      typrisq0,nz0,zi0,nbevt0,idtrunc0,logspecif0,
                      ny0,ns0,nv0,
                      nobs0,nmes0,npm0,nea0,nfix0,refix0,fix0,
                      idmodel0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0,nvalSPLORD0,
                      ind0,nmescur0)
{
  res <- 0
  .Fortran(C_dens,as.double(re0),
           as.double(best0),
           as.double(Y0),as.double(X0),as.double(Tentr0),as.double(Tevt0),as.integer(Devt0),
           as.integer(ideafp0),as.integer(idgfp0),as.integer(idea0),as.integer(idg0),as.integer(idginterac0),as.integer(idgs0),as.integer(idsurv0),as.integer(idsurvs0),as.integer(idsurvsfp0),as.integer(sharedtype0),as.integer(nGK0),
           as.integer(typrisq0),as.integer(nz0),as.double(zi0),as.integer(nbevt0),as.integer(idtrunc0),as.integer(logspecif0),
           as.integer(ny0),as.integer(ns0),as.integer(nv0),
           as.integer(nobs0),as.integer(nmes0),as.integer(npm0),as.integer(nea0),as.integer(nfix0),as.double(refix0),as.integer(fix0),
           as.integer(idmodel0),as.integer(idlink0),as.integer(nbzitr0),as.double(zitr0),as.double(uniqueY0),as.integer(indiceY0),as.integer(nvalSPLORD0),
           as.integer(ind0),as.integer(nmescur0),dens_res=as.double(res))$dens_res
}
