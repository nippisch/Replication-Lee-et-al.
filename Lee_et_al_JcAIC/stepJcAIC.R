#setwd("C:/Users/lynye/Desktop/FU/Diss/EBPModSel/JacobianAIC/LME")
library(MASS)
library(magic)
library(nlme)
#library(mnormt)

JcAIC <- function(lmefit, domain, trans = "no", lambdaopt = NA, m = NA, 
                    y_col = NULL, interval = c(-2, 2), B = NULL){

  if (is.element(trans, c("box.cox", "dual.power", "log.shift"))){
    
    lmefit0 <- lmefit
    fixed <- formula(lmefit0$call$fixed)
    fixedy <- as.formula(update.formula(fixed, as.formula(paste(y_col, "~ ."))))
    
    if (is.na(lambdaopt)){
      lambdaopt <- optimal_parameter(generic_opt,
                                     fixed = fixedy,
                                     smp_data = lmefit0$data,
                                     smp_domains = domain,
                                     transformation = trans,
                                     interval = interval)
    }
    
    cat(paste0("\nOptimal lambda: ", lambdaopt, ", formula = ", fixedy))
    
    y <- lmefit0$data[[y_col]]
    
    
    if (trans == "box.cox") {
      y_trans <- box_cox(lambdaopt, y)$y
      m <- box_cox(lambdaopt, y)$m
    }
    if (trans == "dual.power"){
      y_trans <- dual_power(lambdaopt, y)$y
      m <- dual_power(lambdaopt, y)$m
    }
    if (trans == "log.shift") {
      y_trans <- log_shift(lambdaopt, y)$y
      lambdaopt <- lambdaopt + log_shift(lambdaopt, y)$m
    }
    
    dat <- lmefit0$data
    dat[[ fixed[[2]] ]] <- y_trans
    lmefit <- update(lmefit0, data = dat) 
  }
  
  X <- model.matrix(as.formula(lmefit$call$fixed), data = lmefit$data)
  unique(lmefit$data[[domain]])
  dom <- unique(lmefit$data[[domain]])
  D <- length(unique(lmefit$data[[domain]])) # Number of domains in data
  N <- nobs(lmefit)
  p <- ncol(X)
  
  Z <- matrix(rep(1,
                  length(lmefit$data[[domain]][lmefit$data[[domain]] == dom[1]])),
              ncol = 1) # start Z
  
  for (i in 2:D) {
    Znext <- matrix(rep(1,
                        length(lmefit$data[[domain]][lmefit$data[[domain]] == dom[i]])),
                    ncol = 1)
    Z <- adiag(Z, Znext)
  }
  
  # variance of random effect, u
  sigma_u2 <- as.numeric(nlme::VarCorr(lmefit)[1,1])
  # variance of individual error, e
  sigma_e2 <- lmefit$sigma^2
  
  y <- lmefit$data[[lmefit$call$fixed[[2]]]]
  
  mu <- fitted(lmefit)

# Since conV_y is a diagonal matrix, the inverse of conV_y can be defined as follows:
  inv_conV_y <- diag((sigma_e2)^(-1), ncol = N, nrow = N) 

# Conditional loglikelihood  
  clogLik <- -0.5*(N*log(2*pi) + N*log(sigma_e2) + t(y - mu) %*%
          inv_conV_y %*% (y - mu))

# cAIC Vaida & Blanchard 2005 and conventional cAIC by Greven & Kneib 2010-----
  if (is.null(B)){
    # variance matrix of random effect u, Sigma_u
    Varu_zero <- sigma_u2 * diag(x = 1, nrow = D, ncol = D, names = TRUE) # D_0 in in Vaida & Blanchard 2005
    Varu_scaled <- Varu_zero/sigma_e2 # Dhat_star in Vaida & Blanchard 2005
    
    # in Vaida & Blanchard 2005
    # rho = tr(H), H = (X Z) (M'M)^-1 (X Z)' with M'M = (X'X  X'Z / Z'X  Z'Z + Varu_scaled)
    
    X_Z <- cbind(X, Z) # block matrix (X Z)
    
    XTX <- t(X) %*% X
    ZTX <- t(Z) %*% X
    XTZ <- t(X) %*% Z
    ZTZ_Varu_scaled <- t(Z) %*% Z + solve(Varu_scaled)
    
    dim(XTX)
    dim(XTZ)
    XTX_XTZ <- cbind(XTX, XTZ)
    dim(XTX_XTZ)
    
    dim(ZTX)
    dim(ZTZ_Varu_scaled)
    ZTX_ZTZ_Varu_scaled <- cbind(ZTX, ZTZ_Varu_scaled)
    dim(ZTX_ZTZ_Varu_scaled)
    
    MTM <- rbind(XTX_XTZ, ZTX_ZTZ_Varu_scaled)
    dim(MTM)
    #start_time <- Sys.time()
    H1 <- (X_Z) %*% solve(MTM) %*% t(X_Z)
    #end_time <- Sys.time()
    #time_dur <- end_time - start_time
    #time_dur
    
    
    rhohat <- sum(diag(H1)) # equivalent to conventional cAIC by Greven & Kneib 2010
    
    # Penalty term of cAIC with ML with unknown sigma_e
    K <- (N * (N - p - 1) * (rhohat + 1)) / ((N - p) * (N - p - 2)) +
      (N * (p + 1)) / ((N - p) * (N - p - 2))
    
    cAIChat <- -2 * clogLik + 2 * BC
 # No trafo (JcAIC = cAIC)  
    if (is.null(trans)){JcAIC <- cAIChat}
 # log trafo (JcAIC without bootstrap bias correction)
    if (trans == "log") {
      if(min(lmefit$data[[y_col]]) <= 0) {m <- abs(min(lmefit$data[[y_col]])) + 1}
      else {m = 0}
      logJacobian <- - sum(log(lmefit$data[[y_col]] + m))
      JcAIC <- cAIChat - 2 * logJacobian
    }
 # box.cox trafo (JcAIC without bootstrap bias correction)
   if (trans == "box.cox"){
     logJacobian <- (lambdaopt - 1) * sum(log(lmefit$data[[y_col]] + m))
     JcAIC <- cAIChat - 2 * logJacobian
   }
 # dual.power trafo (JcAIC without bootstrap bias correction)
   if (trans == "dual.power"){
     logJacobian <- sum(log((lmefit$data[[y_col]] + m)^(lambdaopt - 1) +
                              (lmefit$data[[y_col]] + m)^(-lambdaopt - 1))) - log(2)
     JcAIC <- cAIChat - 2 * logJacobian
   } 
 # log.shift trafo (JcAIC without bootstrap bias correction)
   if (trans == "log.shift"){
     logJacobian <-  - sum(log(lmefit$data[[y_col]])) #+ 1 + abs(min(0, lmefit$data[[y_col]])) )
     JcAIC <- cAIChat - 2 * logJacobian
    }   
  } else{
#-------------------------------------------------------------------------------
# JcAIC with bootstrap bias correction  
    bootBC <- NA
    logJacobian <- NA
# No trafo (JcAIC = cAIC)
    if (trans == "no"){
      #bootBC <- bootbiascorrect_cond_normal(lmefit, B = B)$B1
      bootBC <- boot_biascorrection_DD_normal(lmefit, Z = Z, X = X, D = D, B = B)$B1
      cAIChat_B <- -2 * clogLik + 2 * bootBC
      JcAIC <- cAIChat_B # added by Yeonjoo
    } 
# log trafo (JcAIC with bootstrap bias correction)
   if (trans == "log"){
        if(min(lmefit$data[[y_col]]) <= 0) {m <- abs(min(lmefit$data[[y_col]])) + 1}
        else {m = 0}
        logJacobian <- - sum(log(lmefit$data[[y_col]] + m))
        bootBC <- boot_biascorrection_DD_log(lmefit, Z = Z, X = X, D = D, m = m, 
                                               B = B, y_col)$B1
        JcAIC_B <- -2 * clogLik - 2 * logJacobian + 2 * bootBC
        JcAIC <- JcAIC_B # added by Yeonjoo
      }
    if (trans == "box.cox"){
        logJacobian <- (lambdaopt - 1) * sum(log(lmefit$data[[y_col]] + m))
        bootBC <- boot_biascorrection_DD(lmefit, Z = Z, X = X, D = D, 
                                         lambdaopt = lambdaopt, m = m, 
                                         interval = interval, B = B, y_col = y_col)$B1
        JcAIC_B <- -2 * clogLik - 2 * logJacobian + 2 * bootBC
        JcAIC <- JcAIC_B # added by Yeonjoo
    }
  }

  
  cAIC_out <- list(clogLik = as.numeric(clogLik), 
                   bootBC = as.numeric(bootBC),
                   JcAIC = as.numeric(JcAIC),
                   trans = trans,
                   logJacobian = as.numeric(logJacobian),
                   lambdaopt = lambdaopt,
                   lmefit = lmefit) 
  
  #cat(paste0("LogLikelihood: ", clogLik, ", Jacoabian: ", logJacobian, ", BC: ", bootBC, "\n"))
}

################################################################################
dropterm_JcAIC <- function(object, scope, trace = TRUE, trans = NULL, 
                          interval =  c(-2, 2), B, m = NA, y_col = y_col, ...)
{
  tl <- attr(terms(object$call$fixed), "term.labels") 
  if(missing(scope)) 
    scope <- drop.scope(object$fixed)
  else {
    if(!is.character(scope))
      scope <- attr(terms(update.formula(object$fixed, scope)), "term.labels")
    if(!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L,
                dimnames = list(c("<none>", scope), c("edf", "JcAIC")))
  
  ans[1, ] <- c(extractAIC(object)[1], 
                JcAIC(object, domain = names(object$coefficients$random),
                        trans = trans, interval = interval, B = B, m = m, 
                        y_col = y_col)$JcAIC) ##### ACHTUNG
  
  n0 <- nobs(object)
  
  for(i in seq_len(ns)) {
    tt <- scope[i]
    nfit <- object
    nfit$call$fixed <- update(object$call$fixed, as.formula(paste("~ . -", tt)),
                              evaluate = FALSE) 
    nfit$call$formula <- NULL
    catmessage <- capture.output(nfit <- eval(nfit$call))
    
    ans[i+1, ] <- c(extractAIC(nfit)[1], 
                    JcAIC(nfit, domain = names(nfit$coefficients$random),
                            trans = trans, interval = interval, B = B, m = m, 
                            y_col = y_col)$JcAIC)
    nnew <- nobs(nfit)
    if(all(is.finite(c(n0, nnew))) && nnew != n0)
      stop("number of rows in use has changed: remove missing values?")
  }
  dfs <- ans[1L , 1L] - ans[, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, JcAIC = ans[,2])
  head <- c("Single term deletions", "\nModel:", deparse(formula(object$call$fixed)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

################################################################################
addterm_JcAIC <- function(object, scope, trace = TRUE, trans = NULL, 
                         interval =  c(-2, 2), B, m = NA, y_col = y_col, ...)
{
  if(missing(scope) || is.null(scope)) stop("no terms in scope")
  if(!is.character(scope))
    scope <- add.scope(object$fixed, update.formula(object$fixed, scope))
  
  if(!length(scope))
    stop("no terms in scope for adding to object")
  
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L,
                dimnames = list(c("<none>", scope), c("edf", "JcAIC")))
  ans[1, ] <- c(extractAIC(object)[1], 
                JcAIC(object, domain = names(object$coefficients$random),
                        trans = trans, interval = interval, B = B, m = m, 
                        y_col = y_col)$JcAIC)
  n0 <- nobs(object)
  for(i in seq_len(ns)) {
    tt <- scope[i]
    nfit <- object
    nfit$call$fixed <- update(object$call$fixed, as.formula(paste("~ . +", tt)),
                              evaluate = FALSE)
    nfit$call$formula <- NULL
    catmessage <- capture.output(nfit <- eval(nfit$call))
    
    ans[i+1L, ] <- c(extractAIC(nfit)[1], 
                     JcAIC(nfit, domain = names(nfit$coefficients$random),
                             trans = trans, interval = interval, B = B, m = m, 
                             y_col = y_col)$JcAIC) 
    nnew <- nobs(nfit)
    if(all(is.finite(c(n0, nnew))) && nnew != n0)
      stop("number of rows in use has changed: remove missing values?")
  }
  
  dfs <- ans[, 1L] - ans[1L, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, JcAIC = ans[, 2L])
  head <- c("Single term additions", "\nModel:", deparse(formula(object$call$fixed)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

################################################################################

stepJcAIC <- function(object, scope = NULL, m = NA, y_col,
                      direction = "both", trace = TRUE,
                      steps = 1000, trans = "no", interval = c(-2, 2), B = NULL){
  
  # if ((criteria == "AICc" || criteria == "AICb1" || 
  #     criteria == "AICb2"|| criteria == "KICc" || 
  #    criteria == "KICb1"|| criteria == "KICb2") && 
  #  (is.null(object$model$seed)))  {
  #object$call$seed <- 123
  #catmessage <- capture.output(object <- eval(object$call))
  #cat("Seed in fh object not defined, 123 used as default seed. \n")
  #}

  startobject <- object
  
  cut.string <- function(string) {
    if (length(string) > 1L) 
      string[-1L] <- paste0("\n", string[-1L])
    string
  }
  
  step.results <- function(models, fit, object, catmessage) { 
    change <- vapply(models, "[[" , "change", FUN.VALUE = character(1))
    rdf <- vapply(models, "[[", "df.resid", FUN.VALUE = numeric(1))
    ddf <- c(NA, diff(rdf))
    infcriteria <- NULL
    infcriteria <- vapply(models, "[[", "JcAIC", FUN.VALUE = numeric(1))
    lambdaopt <- NULL
    lambdaopt <- vapply(models, "[[", "lambdaopt", FUN.VALUE = numeric(1))
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
                 "\nInitial Model:", deparse(object$call$fixed), "\nFinal Model:", 
                 deparse(fit$call$fixed), "\n")
    aod <- data.frame(Step = I(change), Df = ddf, criteria = infcriteria, lambda = lambdaopt,
                      check.names = FALSE)  
    attr(aod, "heading") <- heading
    fit$anova <- aod
    #fit$call$MSE <- startobject$call$MSE
    #fit$call$formula <- NULL
    #  invisible(fit <- eval(fit$call))
    #  if (catmessage == TRUE){
    #    cat("\n")
    #    cat("Please note that the model selection criteria are only computed based on 
    #     the in-sample domains. \n \n ")
    #  }
    #  class(fit) <- c("step_fh", "fh", "emdi")
    list(fit = fit, models = models)
  }
  
  Terms <- terms(object)
  object$call$fixed <- object$formula <- Terms
  md <- missing(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  if (missing(scope)) {
    fdrop <- numeric()
    fadd <- attr(Terms, "factors")
    if (md) 
      forward <- FALSE
  }
  else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) 
        attr(terms(update.formula(object$call$fixed, fdrop)), "factors")
      else numeric()
      fadd <- if (!is.null(fadd <- scope$upper)) 
        attr(terms(update.formula(object$call$fixed, fadd)), "factors")
    }
    
  }
  models <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object 
  
  domain <- names(fit$coefficients$random)
  JcAIC_result <- JcAIC(lmefit = fit, domain = domain, trans = trans, m = m, 
                           interval = interval, B = B, y_col = y_col)
  JcAIC <- JcAIC_result$JcAIC
  lambdaopt <- JcAIC_result$lambdaopt
  
  
  edf <- extractAIC(fit)[1]
  
  if (is.na(JcAIC)) 
    stop("AIC is not defined for this model, so 'stepAIC' cannot proceed")
  if (JcAIC == -Inf) 
    stop("AIC is -infinity for this model, so 'stepAIC' cannot proceed")
  
  nm <- 1
  
  if (trace == TRUE) {
    cat("Start: JcAIC", " = ", format(round(JcAIC, 2)), "\n", 
        cut.string(deparse(fit$fixed)), "\n\n", sep = "", "\nlambda: ", lambdaopt)
    flush.console()
  }
  models[[nm]] <- list(df.resid = n - edf, 
                       change = "", 
                       JcAIC = JcAIC, 
                       lambdaopt = lambdaopt) 
  
  while (steps > 0) {
    steps <- steps - 1
    infcriteria <- JcAIC
    ffac <- attr(Terms, "factors")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (backward && length(scope$drop)) {
      aod <- dropterm_JcAIC(object = fit, scope = scope$drop, scale = 0, trans = trans, m = m, 
                           interval = interval, B = B, y_col = y_col)
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L], 
                                        sep = " "))
      
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- addterm_JcAIC(object = fit, scope = scope$add, scale = 0, trans = trans, m = m, 
                             interval = interval, B = B, y_col = y_col)
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], sep = " "))
        aod <- if (is.null(aod)) 
          aodf
        else rbind(aod, aodf[-1, , drop = FALSE])
      }
      attr(aod, "heading") <- NULL
      nzdf <- if (!is.null(aod$Df)) 
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nc <- match("JcAIC", names(aod)) 
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc])
      names(aod) <- c("df", "JcAIC")
      if (trace == TRUE) 
        print(aod[o, ])
      if (o[1L] == 1) 
        break
      change <- rownames(aod)[o[1L]]
    }
    
    fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit <- eval.parent(fit)
    nnew <- nobs(fit, use.fallback = TRUE)
    catmessage <- capture.output(fit <- eval(fit$call))
    
    nnew <- nobs(fit, use.fallback = TRUE) 
    if (all(is.finite(c(n, nnew))) && nnew != n) 
      stop("number of rows in use has changed: remove missing values?")
    Terms <- terms(fit)
    
    JcAIC_result <- JcAIC(lmefit = fit, domain = domain, trans = trans, m = m,
                             interval = interval, B = B, y_col = y_col)
    JcAIC <- JcAIC_result$JcAIC
    lambdaopt <- JcAIC_result$lambdaopt
    
    edf <- extractAIC(fit)[1]
    
    if (trace == TRUE) {
      cat("Start: JcAIC", " = ", format(round(JcAIC, 2)), "\n", 
          cut.string(deparse(fit$fixed)), "\n\n", sep = "", "\nlambda: ", lambdaopt)
      flush.console()
    }
    if (JcAIC >= infcriteria + 1e-07) 
      break
    nm <- nm + 1
    models[[nm]] <- list(df.resid = n - edf, 
                         change = "", 
                         JcAIC = JcAIC, 
                         lambdaopt = lambdaopt,
                         JcAIC_result = JcAIC_result)
  }
  #  if (any(grepl(pattern = c("Please note that the model selection criteria are only computed based on 
  #       the in-sample domains."), x = ""))){
  #    catmessage <- TRUE
  #  } else {catmessage <- FALSE}
  step.results(models = models[seq(nm)], fit, object, catmessage) 
  #list(summary(results), invisible(results))
}

################################################################################
# Test
# load("smpDatLog.RData")
# dat <- Smp[[10]]
# 
# dat$logy <- log(dat$y)
# 
# fitfull <- lme(fixed = logy ~ x1 + x2 + x3 + z1, random = ~ 1|idD, 
#                method = "ML", data = dat)
# 
# library(cAIC4)
# cAIC4::cAIC(fitfull)
# 
# cAIC_VB(lmefit = fitfull, domain = "idD")
# 
# fitnull <- lme(fixed = logy ~ 1, random = ~ 1|idD, 
#                method = "ML", data = dat)
# 
# cAIC4::cAIC(fitnull)
# 
# cAIC_VB(lmefit = fitnull, domain = "idD")
# 
# fitopt <- stepcAIC_VB(fitfull, direction = "both", scope = list(lower = fitnull, upper = fitfull))  
# summary(fitopt)

