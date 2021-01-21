# Create initial model
create_TRMF = function(dataM,weight=1,normalize=c("none","standard","robust","range"),
                       normalize.type = c("global","columnwise","rowwise"),
                       na.action=c("impute","fail")){
  
  dataM = as.matrix(dataM)
  
  # Check NAs
  if(match.arg(na.action) == "fail"){
    if(any(is.na(dataM))){
      stop("Missing values in object")
    }
  }
  
  if(any(is.infinite(dataM))){
    stop("Infinite values in data matrix")
  }
  
  if(any(is.infinite(weight))){
    stop("Infinite values in weights")
  }
  
  # Attributes
  Dims = list(nrows = dim(dataM)[1],ncols = dim(dataM)[2],numTS=0)
  
  # Normalize data
  normalize = match.arg(normalize)
  normalize.type = match.arg(normalize.type)
  NormalizedData = NormalizeMatrix(dataM,method=normalize,type=normalize.type)
  NormalizedData[is.na(dataM)]=0
  
  # Projection for missing values
  HadamardProjection = HadamardProjection4NA(dataM)
  
  # Get data weights
  weight[is.na(weight)]=0
  Weight = weight*HadamardProjection
  
  # Create TRMF object
  trmf_object = list(dataM = dataM,NormalizedData = NormalizedData,HadamardProjection=HadamardProjection,
                     Weight=Weight,Dims = Dims,HasXreg=FALSE)
  
  class(trmf_object) = "TRMF"
  
  return(trmf_object)
  
}

# Add coefficient model
TRMF_coefficients = function(obj,reg_type =c("l2","nnls","constrain","interval","none"),lambda=0.0001){
  
  # check object
  if(is.null(obj)||class(obj) != "TRMF"){
    stop("TRMF_coefficients: Create a valid TRMF object first using TRMF_model()")
  }
  
  if(!is.null(obj$Fm_settings)){
    warning("TRMF_coefficient model already defined, overwriting")
  }
  
  # screen constraint type
  type = match.arg(reg_type)
  if(!(type %in%c("l2","nnls","constrain","interval","none"))){
    stop("TRMF_coefficients: Coefficient regularization type not valid (at least not currently implemented)")
    # could also add interval at some point.
  }
  
  # verify lambda
  if(type=="none"){
    lambda=0
  }
  
  
  # update object
  obj$Fm_Settings = list(type=type,lambda=lambda)
  
  return(obj)
}


# Add slope constraint model
TRMF_trend = function(obj,numTS = 1,order = 1,lambdaD=1,lambdaA=0.0001,weight=1){
  
  # verify object
  if(is.null(obj)||class(obj) != "TRMF"){
    stop("TRMF_trend: Create a valid TRMF object first using TRMF_model()")
  }
  
  if(any(is.infinite(weight))){
    stop("TRMF_trend: Infinite values in weights")
  }
  
  # check inputs
  numTS = as.integer(numTS)
  if(numTS<1){
    return(obj)
  }
  
  
  if((length(lambdaD)!=1)||(length(lambdaA)!=1)){
    stop("TRMF_trend: the regularization parameters (lambda) must be scalars")
  }
  
  # list of rules for Xm
  if(is.null(obj$Xm_models)){
    xm_it = 1
    obj$Xm_models=list()
  }else{
    xm_it = length(obj$Xm_models)+1
  }
  
  
  # Check weight input, make weights have mean 1, put in matrix form
  nrows = obj$Dims$nrows
  if((length(weight)!=1)&&(length(weight)!=nrows)){
    stop("TRMF_trend: weight vector is wrong size")
  }
  
  #sumwt = 1/mean(weight)
  #if(is.infinite(sumwt)||is.na(sumwt)){
  #  stop("TRMF_trend: weight vector not valid")
  #}else{
  #  weight = weight/sumwt
  #}
  
  WeightD = diag(x=lambdaD*weight,nrow=nrows)
  
  # Create finite difference constraint
  WeightD = diag(x=lambdaD*weight,nrow=nrows)
  Dmat = FiniteDiffM(nrows,order)
  Dmat = WeightD%*%Dmat
  
  # Overall regularization
  WeightA = diag(x=lambdaA,nrow=nrows)
  
  # Create Xm object
  XmObj = list(Rm = Dmat,WA = WeightA)
  XmObj$model =list(type = "trend",order=order,numTS=numTS)
  XmObj$model$name = paste("Order_",order," trend with ",numTS," latent time series",sep="",collapse="")
  XmObj$model$colnames = paste("D",round(order,2),"(",1:numTS,")",sep="")
  obj$Xm_models[[xm_it]] = XmObj
  obj$Dims$numTS=obj$Dims$numTS+numTS
  return(obj)
}


# No temporal structure
TRMF_simple = function(obj,numTS = 1,lambdaA=0.0001,weight=1){
  
  # verify object
  if(is.null(obj)||class(obj) != "TRMF"){
    stop("TRMF_trend: Create a valid TRMF object first using TRMF_model()")
  }
  
  if(any(is.infinite(weight))){
    stop("TRMF_simple: Infinite values in weights")
  }
  
  # check inputs
  numTS = as.integer(numTS)
  if(numTS<1){
    return(obj)
  }
  
  
  if(length(lambdaA)!=1){
    stop("TRMF_simple: the regularization parameter (lambda) must be scalar")
  }
  
  # list of rules for Xm
  if(is.null(obj$Xm_models)){
    xm_it = 1
    obj$Xm_models=list()
  }else{
    xm_it = length(obj$Xm_models)+1
  }
  
  
  # Check weight input, make weights have mean 1, put in matrix form
  nrows = obj$Dims$nrows
  if((length(weight)!=1)&&(length(weight)!=nrows)){
    stop("TRMF_simple: weight vector is wrong size")
  }
  
  WeightA = diag(x=lambdaA*weight,nrow=nrows)
  
  
  
  # Overall regularization
  WeightD = diag(x=0,nrow=nrows) # this is here to make consistent with other models
  
  # Create Xm object
  XmObj = list(Rm = WeightD ,WA = WeightA)
  XmObj$model =list(type = "simple",order=0,numTS=numTS)
  XmObj$model$name = paste0("L2 regularized with ",numTS," latent time series")
  XmObj$model$colnames = paste("L2(",1:numTS,")",sep="")
  obj$Xm_models[[xm_it]] = XmObj
  obj$Dims$numTS=obj$Dims$numTS+numTS
  return(obj)
}

# Add seasonal random walk model
TRMF_seasonal = function(obj,numTS = 1,freq = 12,sumFirst=FALSE,lambdaD=1,lambdaA=0.0001,weight=1){
  
  # verify object
  if(is.null(obj)||class(obj) != "TRMF"){
    stop("TRMF_seasonal: Create a valid TRMF object first using TRMF_model()")
  }
  
  if(any(is.infinite(weight))){
    stop("TRMF_seasonal: Infinite values in weights")
  }
  
  # check inputs
  numTS = as.integer(numTS)
  if(numTS<1){
    return(obj)
  }
  
  
  if(round(freq) != freq){
    message("TRMF_seasonal: Non-integer frequencies (freq) currently rounded to nearest integer")
  }
  if(freq<1){
    stop("TRMF_seasonal: freq value not valid")
  }
  if(freq==1){
    message("TRMF_seasonal: lag = freq, consider using TRMF_trend() instead")
  }
  
  if((length(lambdaD)!=1)||(length(lambdaA)!=1)){
    stop("TRMF_seasonal: the regularization parameters (lambda) must be scalars")
  }
  
  # list of rules for Xm
  if(is.null(obj$Xm_models)){
    xm_it = 1
    obj$Xm_models=list()
  }else{
    xm_it = length(obj$Xm_models)+1
  }
  
  # Check weight input, make weights have mean 1, put in matrix form
  nrows = obj$Dims$nrows
  if((length(weight)!=1)&&(length(weight)!=nrows)){
    stop("TRMF_seasonal: weight vector is wrong size")
  }
  
  
  WeightD = diag(x=lambdaD*weight,nrow=nrows)
  
  # Create finite difference constraint
  WeightD = diag(x=lambdaD*weight,nrow=nrows)
  Dmat = Seasonal_DM(nrows,lag=freq,sumFirst=sumFirst)
  Dmat = WeightD%*%Dmat
  
  # Overall regularization
  WeightA = diag(x=lambdaA,nrow=nrows)
  
  # Create Xm object
  XmObj = list(Rm = Dmat,WA = WeightA)
  XmObj$model =list(type = "seasonal",freq=freq,numTS=numTS)
  XmObj$model$name = paste("Frequency = ",freq," seasonal random walk with ",numTS," latent time series",sep="",collapse="")
  XmObj$model$colnames = paste("L",freq,"(",1:numTS,")",sep="")
  obj$Xm_models[[xm_it]] = XmObj
  obj$Dims$numTS=obj$Dims$numTS+numTS
  return(obj)
}

# Add an Auto-Regressive model
TRMF_ar = function(obj,numTS = 1,AR,lambdaD=1,lambdaA=0.0001,weight=1){
  
  # verify object
  if(is.null(obj)||class(obj) != "TRMF"){
    stop("TRMF_AR: Create a valid TRMF object first using TRMF_model()")
  }
  
  if(any(is.infinite(weight))){
    stop("TRMF_ar: Infinite values in weights")
  }
  
  # check inputs
  numTS = as.integer(numTS)
  if(numTS<1){
    return(obj)
  }
  
  
  
  if((length(lambdaD)!=1)||(length(lambdaA)!=1)){
    stop("TRMF_AR: the regularization parameters (lambda) must be scalars")
  }
  
  # list of rules for Xm
  if(is.null(obj$Xm_models)){
    xm_it = 1
    obj$Xm_models=list()
  }else{
    xm_it = length(obj$Xm_models)+1
  }
  
  # Check weight input, make weights have mean 1, put in matrix form
  nrows = obj$Dims$nrows
  if((length(weight)!=1)&&(length(weight)!=nrows)){
    stop("TRMF_AR: weight vector is wrong size")
  }
  
  
  WeightD = diag(x=lambdaD*weight,nrow=nrows)
  
  # Create finite difference constraint
  WeightD = diag(x=lambdaD*weight,nrow=nrows)
  Amat = ARmat(nrows,AR)
  Amat = WeightD%*%Amat
  
  # Overall regularization
  WeightA = diag(x=lambdaA,nrow=nrows)
  
  
  # Create Xm object
  XmObj = list(Rm = Amat,WA = WeightA)
  XmObj$model =list(type = "auto-regressive",parms = AR,numTS=numTS)
  XmObj$model$name = paste("Auto-regressive model of order ",length(AR)," with ",numTS," latent time series",sep="",collapse="")
  XmObj$model$colnames = paste("AR",length(AR),"(",1:numTS,")",sep="")
  obj$Xm_models[[xm_it]] = XmObj
  obj$Dims$numTS=obj$Dims$numTS+numTS
  return(obj)
}

# Add an exponential smoothing model
TRMF_es = function(obj,numTS = 1,alpha=1,es_type=c("single","double"),lambdaD=1,lambdaA=0.0001 ,weight=1){
  
  # verify object
  if(is.null(obj)||class(obj) != "TRMF"){
    stop("TRMF_ES: Create a valid TRMF object first using TRMF_model()")
  }
  
  if(any(is.infinite(weight))){
    stop("TRMF_es: Infinite values in weights")
  }
  
  # check inputs
  numTS = as.integer(numTS)
  if(numTS<1){
    return(obj)
  }
  
  
  if((length(lambdaD)!=1)||(length(lambdaA)!=1)){
    stop("TRMF_ES: the regularization parameters (lambda) must be scalars")
  }
  
  es_type = match.arg(es_type)
  if(es_type=="single"){
    order = 1
    prefix =""
  }else if(es_type=="double"){
    order = 2
    prefix ="D"
  }else{
    stop("TRMF_ES: exponential smoothing type not valid")
  }
  
  # list of rules for Xm
  if(is.null(obj$Xm_models)){
    xm_it = 1
    obj$Xm_models=list()
  }else{
    xm_it = length(obj$Xm_models)+1
  }
  
  
  # Check weight input, make weights have mean 1, put in matrix form
  nrows = obj$Dims$nrows
  if((length(weight)!=1)&&(length(weight)!=nrows)){
    stop("TRMF_ES: weight vector is wrong size")
  }
  
  WeightD = diag(x=lambdaD*weight,nrow=nrows)
  
  # Create finite difference constraint
  WeightD = diag(x=lambdaD*weight,nrow=nrows)
  Dmat = ExpSmMat(nrows,alpha,order)
  Dmat = WeightD%*%Dmat
  
  # Overall regularization
  WeightA = diag(x=lambdaA,nrow=nrows)
  
  # Create Xm object
  XmObj = list(Rm = Dmat,WA = WeightA)
  XmObj$model =list(type = "ES",alpha=alpha,order=order,numTS=numTS)
  XmObj$model$name = paste(es_type," exponential smoothing with ",numTS," latent time series",sep="",collapse="")
  XmObj$model$colnames = paste(prefix,"ES",round(alpha,2),"(",1:numTS,")",sep="")
  obj$Xm_models[[xm_it]] = XmObj
  obj$Dims$numTS=obj$Dims$numTS+numTS
  return(obj)
}

# Add a Regression model
TRMF_regression = function(obj,Xreg,type=c("global","columnwise")){
  
  # verify object
  if(is.null(obj)||class(obj) != "TRMF"){
    stop("TRMF_Regression: Create a valid TRMF object first using TRMF_model()")
  }
  
  if(any(is.infinite(Xreg))){
    stop("TRMF_Regression: infinite values in external regressors")
  }
  
  if(any(is.na(Xreg))){
    stop("TRMF_Regression: Missing values not allowed in external regressors")
  }
  
  nrows = obj$Dims$nrows
  type = match.arg(type)
  
  dimX = dim(Xreg)
  
  # Global case, all the data time series regress off of the same external regressors
  if(type=="global"){
    
    # check dimensions
    if(is.null(dimX)){
      if(length(Xreg) != nrows){
        stop("TRMF_Regression: Xreg dimensions are incompatible with the data")
      }else{
        Xreg = matrix(Xreg,nrow=nrows)
        dimX = dim(Xreg)
      }
      
    }else{
      if(dimX[1] != nrows){
        stop("TRMF_Regression: Xreg dimensions are incompatible with the data")
      }
    }
    if(length(dimX)>2){
      stop("TRMF_Regression: Xreg has more then 2 dimensions, perhaps you meant to use type='columnwise'")
    }
    
    # only allow 1 global regression model
    if(!is.null(obj$GlobalXReg)){
      warning("TRMF_Regression: A global external regressor model has already been defined, over-writing...")
    }
    
    # store
    obj$xReg_model$GlobalXReg = Xreg
    
    # create name
    cnames = colnames(Xreg)
    if(is.null(cnames)){
      cnames = paste("gXreg(",1:dimX[2],")",sep="")
    }else{
      cnames = paste("gXreg(",cnames,")",sep="")
    }
    obj$xReg_model$gxname = cnames
    
    
  }else if(type=="columnwise"){   # local case, each data time series gets its own regressor
    
    
    # check dimensions
    if(is.null(dimX)){
      stop("TRMF_Regression: Xreg is not matrix or array, perhaps you meant to use type='global'")
    }
    if(dimX[1] != nrows){
      stop("TRMF_Regression: number of rows of Xreg do not match number of rows of data")
    }
    if(dimX[2] != obj$Dims$ncols){
      stop("TRMF_Regression: number of columns of Xreg do not match number of columns of data")
    }
    
    # only allow 1 column-wise regression model (can implement more later)
    if(!is.null(obj$ColumnwiseXReg)){
      warning("TRMF_Regression: A columnwise external regressor model has already been defined, over-writing...")
    }
    
    # store as array
    if(is.na(dimX[3])){dimX[3]=1}
    obj$xReg_model$ColumnwiseXReg = array(Xreg,dimX)
    
    # create names
    obj$xReg_model$cxname = paste("cXreg(",1:dimX[3],")",sep="")
    
    # Old attempt, not really possible to get names to line up in this case
    #cnames = dimnames(Xreg)[[2]]
    #if(is.null(cnames)){
    #  cnames = paste("cXreg(",1:dimX[3],sep=")")
    #}else{
    #  cnames = paste("cXreg(",cnames,sep=")")
    #}
    #obj$xReg_model$cxname = cnames
    
    
  }
  
  return(obj)
}

# Fit TRMF model
#train_TRMF = function(obj,numit=10){
  
train.TRMF = function(x,numit=10,...){
  
  # check object to see if it valid object
  if(is.null(obj)||class(obj) != "TRMF"){
    # should never get here...
    stop("Not a valid TRMF object")
  }
  
  
  # if no models have been added, create a simple one
  if(is.null(obj$Fm_Settings)){
    obj = TRMF_coefficients(obj)
  }
  
  if(is.null(obj$Xm_models)){
    obj = TRMF_simple(obj)
  }
  
  # Set up stuff
  localEnv = list2env(obj,envir = new.env())
  Create_Xreg_Stuff(localEnv)
  Create_Fm_Stuff(localEnv)
  Create_Xm_Stuff(localEnv)
  
  # Run ALS iteration
  InitializeALS(localEnv)
  for(k in 1:numit){
    Get_XReg_fit(localEnv)
    FitXm(localEnv)
    FitFm(localEnv)
  }
  
  
  # get fit
  FitAll(localEnv)
  
  # format back as object
  newobj=as.list(localEnv)
  class(newobj) = class(obj)
  return(newobj)
  
}

impute_TRMF= function(obj){
  
  # check object to see if it valid object
  if(is.null(obj)||class(obj) != "TRMF"){
    stop("Not a valid TRMF object")
  }
  
  # Make sure everything we need is there
  if(is.null(obj$Fit)){
    stop("Train TRMF model first using 'train_TRMF()'")
  }
  
  # find missing values
  newM = obj$dataM
  ind = which(is.na(newM))
  newM[ind] = obj$Fit$fitted[ind]
  
  return(newM)
  
}
