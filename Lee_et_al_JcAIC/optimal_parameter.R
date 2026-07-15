
optimal_parameter <- function(generic_opt,
                              fixed,
                              smp_data,
                              smp_domains,
                              transformation,
                              interval = c(-3,3)) {
  
  if(transformation != "no" &&
     transformation != "log") {
    # no lambda -> no estimation -> no optmimization
    
    # Estimation of optimal lambda parameters
    optimal_parameter <- optimize(generic_opt,
                                  fixed          = fixed,
                                  smp_data       = smp_data,
                                  smp_domains    = smp_domains,
                                  transformation = transformation,
                                  interval       = interval,
                                  maximum        = FALSE
    )$minimum
    
  } else {
    optimal_parameter <- NULL
  }
  
  return(optimal_parameter)
} # End optimal parameter


# Internal documentation -------------------------------------------------------

# Function generic_opt provides estimation method reml to specifiy
# the optimal parameter lambda. Here its important that lambda is the
# first argument because generic_opt is given to optimize. Otherwise,
# lambda is missing without default.

generic_opt <- function(lambda,
                        fixed,
                        smp_data,
                        smp_domains,
                        transformation
) {
  
  
  #Definition of optimization function for finding the optimal lambda
  #Preperation to easily implement further methods here
  optimization <- if (TRUE) {
    reml(fixed          = fixed,
         smp_data       = smp_data,
         smp_domains    = smp_domains,
         transformation = transformation,
         lambda         = lambda
    )
  }
  return(optimization)
}



# REML method ------------------------------------------------------------------

reml <- function(fixed          = fixed,
                 smp_data       = smp_data,
                 smp_domains    = smp_domains,
                 transformation = transformation,
                 lambda         = lambda
) {
  
  sd_transformed_data <- std_data_transformation(fixed          = fixed,
                                                 smp_data       = smp_data,
                                                 transformation = transformation,
                                                 lambda         = lambda
  )
  model_REML <- NULL
  try(model_REML <- nlme::lme(fixed     = fixed,
                        data      = sd_transformed_data,
                        random    = as.formula(paste0("~ 1 | as.factor(", smp_domains, ")")),
                        method    = "REML",
                        keep.data = FALSE), silent = TRUE)
  if(is.null(model_REML)) {
    try(model_REML <- nlme::lme(fixed     = fixed,
                          data      = sd_transformed_data,
                          random    = as.formula(paste0("~ 1 | as.factor(", smp_domains, ")")),
                          method    = "REML",
                          control = nlme::lmeControl(opt = "optim"),
                          keep.data = FALSE), silent = TRUE)
    if(is.null(model_REML)) {
      stop("The likelihood does not converge. One reason could be that the
         interval for the estimation of an optimal transformation parameter is
         not appropriate. Try another interval. See also help(ebp).")
    } else {
      model_REML <- model_REML
    }
  } 
  
  log_likelihood <- -logLik(model_REML)
  return(log_likelihood)
}


std_data_transformation <- function(fixed=fixed,
                                    smp_data,
                                    transformation,
                                    lambda) {
  
  y_vector <- as.matrix(smp_data[paste(fixed[2])])
  
  std_transformed <- if (transformation == "box.cox"){
    as.data.frame(box_cox_std(y = y_vector, lambda = lambda))
  } else if (transformation == "dual.power") {
    as.data.frame(dual_power_std(y = y_vector, lambda = lambda))
  } else if (transformation == "log.shift") {
    as.data.frame(log_shift_std(y = y_vector, lambda = lambda))
  } else if (transformation == "log") {
    smp_data[paste(fixed[2])]
  } else if (transformation == "no") {
    smp_data[paste(fixed[2])]
  }
  
  smp_data[paste(fixed[2])] <- std_transformed
  return(transformed_data = smp_data)
} # End std_data_transformation

# Standardized transformation: Box Cox

geometric.mean <- function(x) { #for RMLE in the parameter estimation
  
  exp(mean(log(x)))
}

box_cox_std <- function(y, lambda) {
  min <- min(y)
  if (min <= 0) {
    y <- y - min + 1
  }
  
  gm <- geometric.mean(y)
  y <- if (abs(lambda) > 1e-12) {
    y <- (y^lambda - 1) / (lambda * ((gm)^(lambda - 1)))
  } else {
    y <- gm * log(y)
  }
  return(y)
}

dual_power_std <-  function(y, lambda) {
  min <- min(y)
  if (min <= 0) {
    y <- y - min + 1
  }
  
  y <- if (abs(lambda) > 1e-12) {
    geo <- geometric.mean(y^(lambda-1)+y^(-lambda-1))
    y <- ((y^lambda -y^(-lambda))/(2*lambda))*2/geo
  } else {
    y <- geometric.mean(y) * log(y)
  }
  return(y)
}


dual_power <- function(l, y, m = 0){
  m <- 0
  if(min(y) <= 0)
  {
    m = abs(min(y)) + 1
    y = y + m
  }
  if(abs(l)<= 1e-12)
  {
    y <- log(y)
  }
  else
  {
    y <- (y^l - y^(-l)) / (2 * l)
  }
  return(list(y=y, m=m))
}


box_cox <- function(l, y, m=NULL) #Box-Cox transformation (lambda=l)
{
  if(is.null(m))
  {
    m = 0
  }
  if((s=min(y))<=0) #with shift(=m) parameter for making data positive (>0)
  {
    s = abs(s)+1
  }
  else
  {
    s=0
  }
  m=m+s
  
  if(abs(l)<=1e-12) #case lambda=0
  {
    y = log(y+m)
  }
  else
  {
    y = ((y+m)^l-1)/l
  }
  return(list(y = y, m = m)) #return of transformed data and shift (overwriten y)
}

# Back transformation: Box Cox
box_cox_back <- function(y, lambda, shift = 0) {
  
  lambda_cases_back <- function(y, lambda = lambda, shift){
    if (abs(lambda) <= 1e-12) {   #case lambda=0
      y <-  exp(y) - shift
    } else {
      y <- (lambda * y + 1)^(1 / lambda) - shift
    }
    return(y = y)
  }
  y <- lambda_cases_back(y = y, lambda = lambda, shift = shift)
  
  return(y = y)
} #  End box_cox_back


log_shift <- function(l, y, m=NULL) #Log-shift transformation
{
  
  y <- y + l
  if(is.null(m))
  {
    m = 0
  }
  if((s=min(y))<=0) #with shift(=m) parameter for making data positive (>0)
  {
    s = abs(s) + 1
  }
  else
  {
    s=0
  }
  m = m + s
  y = log(y + m)
  
  return(list(y = y, m = m))
}

log_shift_std <- function(y, lambda, m) # standardized Log-shift
{
  y <- y + lambda
  
  min <- min(y)
  if (min <= 0) {
    y <- y - min + 1
  }
  
  gm <-  geometric.mean(y)
  y <- gm * log(y)
}

bound_estimation = function(y , max_range_bounds = NULL, m=NULL) # bounds for the estimation parameter
{
  
  #for log shift enshure that the data are positive. The other transformations use the given interval
  span=range(y)
  if( (span[1]+1) <= 1)
  {
    lower = abs(span[1])+1
  }
  else
  {
    lower = 0
  }
  upper = diff(span)
  
  return(c(lower,upper))
  
}
