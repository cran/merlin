#' merlin - Mixed Effects Regression for Linear, Nonlinear and User-defined models
#'
#' merlin fits linear, non-linear and user-defined mixed effects regression models. merlin can fit
#' multivariate outcome models of any type, each of which could be repeatedly measured
#' (longitudinal), with any number of levels, and with any number of random effects at each level.
#' Standard distributions/models available include the Bernoulli, Gaussian, Poisson, beta,
#' negative-binomial, and time-to-event/survival models include the exponential, Gompertz,
#' Weibull, Royston-Parmar, and general hazard model. merlin provides a
#' flexible predictor syntax, allowing the user to define variables, random effects, spline and
#' fractional polynomial functions, functions of other outcome models, and any interaction
#' between each of them. Non-linear and time-dependent effects are seamlessly incorporated into
#' the predictor. merlin allows multivariate normal random effects, which are integrated out
#' using Gaussian quadrature or Monte-Carlo integration. Relative survival (excess hazard) models
#' are supported. Utility functions are provided to allow user-defined models to be specified,
#' in conjunction with the complex predictor.
#'
#' @author Emma C. Martin, Alessandro Gasparini and Michael J. Crowther
#'
#' @param model specify the fixed and random elements for each model outcome.
#'   Where there are multiple outcomes, the models should be specified in a list.
#'   Each model should be specified as a formula (e.g. \code{y~x}). A number of different element types
#'   can be specified, including
#'   \itemize{
#'     \item varname - an independent variable from the data set
#'     \item random-effects - a random-effect at the cluster level can be specified
#'       using \code{M#[cluster level]}, for example \code{M1[id]} would define a random intercept
#'       at the ID level. Each independent random-effect should be given a unique name, if two
#'       random-effects are given the same name they will be treated as shared random-effects.
#'     \item \code{rcs()} - restricted cubic spline terms, this option can be used to include a
#'       restricted cubic spline function, with the degrees of freedom (number of spline terms)
#'       specified using the \code{df} sub-option, with the boundary knots assumed to be at the
#'       minimum and maximum of the variable, with internal knots placed at equally spaced centiles. Other
#'       default options \code{orthog=TRUE}, which by default orthogonalises the splines,
#'       \code{log=FALSE}, which can be used to calculate splines of the log of the variable
#'       and \code{event=F}, which can be used to calculate the internal knots based only on observations
#'       that experienced the event of interest (for survival models).
#'     \item \code{srcs()} is a shorthand element, equivalent to \code{rcs(...,log=T,event=T)}, for use with
#'     survival models.
#'     \item \code{fp()} - fractional polynomials of order 1 or 2 can be specified, the
#'       sub-option \code{powers} is used to specify the powers of the fp model.
#'     \item \code{EV[depvar]} - the expected value of the response of a submodel
#'     \item \code{dEV[depvar]} - the first derivative with respect to time of the expected value
#'     of the response of a submodel
#'     \item \code{d2EV[depvar]} - the second derivative with respect to time of the expected value
#'     of the response of a submodel
#'     \item \code{iEV[depvar]} - the integral with respect to time of the expected value
#'     of the response of a submodel
#'     \item \code{XB[depvar]} - the expected value of the complex predictor of a submodel
#'     \item \code{dXB[depvar]} - the first derivative with respect to time of the expected value
#'     of the complex predictor of a submodel
#'     \item \code{d2XB[depvar]} - the second derivative with respect to time of the expected value
#'     of the complex predictor of a submodel
#'     \item \code{iXB[depvar]} - the integral with respect to time of the expected value
#'     of the complex predictor of a submodel
#'     \item \code{bhazard(varname)} - invokes a relative survival (excess hazard) model. \code{varname}
#'     specifies the expected hazard rate at the event time.
#'     \item \code{exposure(varname)} - include \code{log(varname)} in the linear predictor, with a
#'     coefficient of 1. For use with \code{family = "poisson"}.
#'     \item \code{ap(#)} - defines the number of ancillary parameters. Used with \code{family="user"}.
#'   }
#' @param family a vector of strings specifying the family for each outcome specified in model.
#' The currently available models include,
#' \itemize{
#'   \item \code{gaussian} - Gaussian distribution
#'   \item \code{bernoulli} - Bernoulli distribution
#'   \item \code{poisson} - Poisson distribution
#'   \item \code{beta} - Beta distribution
#'   \item \code{negbinomial} - Negative binomial distribution
#' }
#' with survival models including,
#' \itemize{
#'   \item \code{exponential} - exponential distribution
#'   \item \code{weibull} - Weibull distribution
#'   \item \code{gompertz} - Gompertz distribution
#'   \item \code{rp} - Royston-Parmar model (complex predictor on the log cumulative hazard scale)
#'   \item \code{loghazard} - general log hazard model (complex predictor on the log hazard scale)
#' }
#' and user-defined,
#' \itemize{
#'   \item \code{user} - fit a user-defined model, which should be written using \code{merlin}'s
#'   utility functions. The name of your user-defined function should be passed through the
#'   \code{userf} option.
#'   \item \code{null} - is a convenience tool to define additional complex predictors, that
#'   do not contribute to the log likelihood. For use with \code{family="user"}.
#' }
#' @param link string vector defining the link functions for each model. Default is \code{"identity"} for
#' all models. If specified, you must define a link function for all submodels. Options include
#' \code{"identity"}, \code{"log"} and \code{"logit"}.
#' @param userf string vector defining the name of the user-written functions for each \code{family="user"}.
#' Each function must be in memory and should return the observation level contribution to the log-likelihood.
#' @param timevar specifies the variable which represents time, this is necessary when a function of time
#' is used in the linear predictor of a survival model as it may interact with other elements of the model.
#' @param intmethod method used for numerically integrating out the random-effects in order to
#' calculate the likelihood for a mixed effects model which includes random effects.
#' Options include \code{ghermite} for non-adaptive Gauss-Hermite quadrature, \code{halton} for
#' Monte-Carlo integration using Halton sequences, or \code{sobol} for Monte-Carlo integration
#' using Sobol sequences, or \code{mc} for standard Monte-Carlo integration using normal draws. The
#' default is \code{ghermite}. Level-specific integration techniques can
#' be specified, for example, with a three level model, we may use Gauss-Hermite quadrature at the
#' highest level, and Monte-Carlo integration with Halton sequences at level 2, using
#' \code{intmethod = c("ghermite","halton")}.
#' @param data a data frame containing all variables required for fitting the model. Can be a tibble object.
#' @param covariance the structure of the variance-covariance matrix can be varied, the default is
#' \code{diagonal} where all diagonal elements of the variance-covariance matrix are estimated uniquely,
#' \code{identity} assumes all diagonal elements are equal and \code{unstructured} estimates all
#' elements of the variance-covariance matrix.
#' @param ip an optional vector of integers specifying the number of integration points
#' to be used when integrating out the random effects. A different number of \code{ip} can be
#' specified for each \code{level} (from highest to lowest level). If only a single number is given then
#' it will assume that number of integration points at all levels. Default is \code{ip = 7} for each
#' level using Gauss-Hermite quadrature, or \code{ip = 100} for each level using Monte-Carlo
#' integration.
#' @param levels if the model contains random-effects then a vector giving the order of levels must
#' be specified, from the highest level to the lowest, e.g. \code{levels=c("practice","id")}.
#' @param from this is an optional argument giving the initial values for the full parameter
#' vector, for more details on how to specify the initial estimates see the vignette.
#' @param sweights Not documented.
#' @param debug Not documented.
#' @param verbose Not documented.
#' @param predict Not documented.
#' @param predtype Not documented.
#' @param predmodel Not documented.
#' @param causes Not documented.
#'
#' @seealso \code{\link{predict.merlin}}
#' @seealso \code{\link{merlin_util_depvar}}, \code{\link{merlin_util_timevar}},
#' \code{\link{merlin_util_xzb}}, \code{\link{merlin_util_xzb_mod}}
#' \code{\link{merlin_util_xzb_deriv}}, \code{\link{merlin_util_xzb_deriv_mod}}
#' \code{\link{merlin_util_xzb_deriv2}}, \code{\link{merlin_util_xzb_deriv2_mod}}
#' \code{\link{merlin_util_xzb_integ}}, \code{\link{merlin_util_xzb_integ_mod}}
#' \code{\link{merlin_util_ev}}, \code{\link{merlin_util_ev_mod}}
#' \code{\link{merlin_util_ev_deriv}}, \code{\link{merlin_util_ev_deriv_mod}}
#' \code{\link{merlin_util_ev_deriv2}}, \code{\link{merlin_util_ev_deriv2_mod}}
#' \code{\link{merlin_util_ev_integ}}, \code{\link{merlin_util_ev_integ_mod}}
#' \code{\link{merlin_util_ap}}, \code{\link{merlin_util_ap_mod}}
#'
#' @references Crowther MJ. Extended multivariate generalised linear and non-linear mixed effects models. \url{https://arxiv.org/abs/1710.02223}
#' @references Crowther MJ. merlin - a unified framework for data analysis and methods development in Stata. \url{https://arxiv.org/abs/1806.01615}
#' @references Martin EC, Gasparini A, Crowther MJ. merlin - an R package for mixed effects regression of linear, non-linear and user-defined models.
#'
#' @examples
#' \donttest{
#' library(merlin)
#' data(pbc.merlin, package = "merlin")
#'
#' # Linear fixed-effects model
#' merlin(logb ~ year,
#'        family = "gaussian",
#'        data = pbc.merlin)
#'
#' # Linear mixed-effects model with random intercept and slope at ID level
#' merlin(logb ~ year + M1[id] * 1 + year:M2[id] * 1,
#'        family = "gaussian",
#'        levels = "id",
#'        data = pbc.merlin)
#'
#' # Joint longitudinal and survival model with shared random effects
#' merlin(model = list(logb ~ year + M1[id] * 1,
#'                     Surv(stime, died) ~ trt + M1[id]),
#'        family = c("gaussian", "weibull"),
#'        levels = "id",
#'        data = pbc.merlin)
#'
#' # Joint longitudinal and survival model with expected value
#' merlin(model = list(logb ~ year + M1[id] * 1,
#'                     Surv(stime, died) ~ trt + EV[logb]),
#'        family = c("gaussian", "weibull"),
#'        levels = "id",
#'        timevar = c("year","stime"),
#'        data = pbc.merlin)
#'
#' # Gaussian distribution - implemented as a user-written family
#' logl_gaussian <- function(gml)
#' {
#'   y    <- merlin_util_depvar(gml)
#'   xzb  <- merlin_util_xzb(gml)
#'   se   <- exp(merlin_util_ap(gml,1))
#'
#'   mu   <- (sweep(xzb,1,y,"-"))^2
#'   logl <- ((-0.5 * log(2*pi) - log(se)) - (mu/(2 * se^2)))
#'   return(logl)
#' }
#'
#' merlin(logb ~ year + ap(1), family = "user", data = pbc.merlin,
#'                                   userf = "logl_gaussian")
#'
#' # 3-level Weibull model
#' merlin(Surv(stime1,dead1) ~ age + M1[id1]*1 + M2[id2]*1,
#'        levels=c("id1","id2"), family="weibull", data=sim3)
#' }

#' @export
merlin <- function(model,
                   from=NULL,
                   ip=NULL,
                   family="gaussian",
                   link=NULL,
                   timevar=NULL,
                   covariance="diagonal",
                   intmethod="ghermite",
                   data,
                   userf=NULL,
                   sweights=NULL,
                   levels=NULL,
                   debug=FALSE,
                   verbose=FALSE,
                   predict=FALSE,
                   predtype=NULL,
                   predmodel=NULL,
                   causes=NULL)
{


    # Arguments checks
    if (is.list(model)==F) model <- list(model)
    merlin_error_check(model=model,data=data,timevar=timevar,
                       family=family,link=link,intmethod=intmethod,
                       covariance=covariance,
                       levels=levels,sweights=sweights)

    # Turn tibbles into plain data.frame objects
    data.name <- substitute(data)
    class(data) <- "data.frame"

    est <- merlinEst(model,from,ip,family,link,timevar,
                      covariance,intmethod,data,userf,sweights,
                      levels,debug,verbose,predict,predtype,
                      predmodel,causes)

    if (predict==FALSE) {
        est$call   <- match.call()
        est$data   <- data.name
        class(est) <- "merlin"
    }
    est
}

# MERLIN
merlinEst <- function(model,
                      from=NULL,
                      ip=7,
                      family="gaussian",
                      link=NULL,
                      timevar=NULL,
                      covariance="diagonal",
                      intmethod="ghermite",
                      data,
                      userf=NULL,
                      sweights=NULL,
                      levels=NULL,
                      debug=FALSE,
                      verbose=FALSE,
                      predict=FALSE,
                      predtype=NULL,
                      predmodel=NULL,
                      causes=NULL)
{

    #setup
    gml <- merlin_setup(model=model,intmethod=intmethod,ip=ip,data=data,timevar=timevar,
                        family=family,link=link,covariance=covariance,userf,
                        from=from,levels=levels,debug=debug,sweights=sweights)

    # initial values
    if (length(from) > 1)  b <- from
    else                   b <- merlin_initial_values(gml,model)

    #fit full model or get predictions
    gml$predict <- predict

    if (predict==FALSE) {

        result <- stats::optim(par=b,merlin_gf,gml=gml,hessian=T,control = list(maxit=5000))
        #gr=merlin_gf_deriv
        if (debug==TRUE) {
            print("optimisation complete")
            print(result$par)
            print(result$hessian)
            print(result$convergence)
        }

    }
    else {
        gml$par      = from
        gml$modelind = gml$modtouse = predmodel
        if (length(causes)) gml$causes = causes
        pred = merlin_predict(gml,predtype)
        return(pred)
    }

    # create coeff table
    btab           = matrix(NA,ncol=6,nrow=length(result$par))
    colnames(btab) = c("Coef.","Std. Error","z","P>|z|","[95% Conf.","Interval]")
    rownames(btab) = gml$labelso

    btab[,1] = result$par
    btab[,2] = round(sqrt(diag(solve(result$hessian))),digits = 7)
    btab[,3] = round(btab[,1] / btab[,2], digits = 2)
    btab[,5] = round(btab[,1] - stats::qnorm(0.975) * btab[,2],digits = 7)
    btab[,6] = round(btab[,1] + stats::qnorm(0.975) * btab[,2],digits = 7)
    btab[,4] = round(2*stats::pnorm(-abs(btab[,3])),digits = 3)

    list( "convergence"    = result$convergence,
          "data"           = data,
          "family"         = family,
          "link"           = link,
          "loglikelihood"  = - result$value,
          "levels"         = gml$levels,
          "responsevar"    = gml$y,
          "gamma"          = gamma,
          "coefftable"     = btab,
          "intmethod"      = intmethod,
          "ip"             = gml$ip,
          "par"            = result$par,
          "Nmodels"        = gml$Nmodels
          )
}

#' @export
print.merlin <- function(x, ...) {

    cat("Merlin: mixed-effects model \n")
    cat("Data: ")
    print(x$data)
    cat("Log likelihood = ")
    cat(x$loglikelihood)
    cat("\n")
    cat("\n")
    cat("Results:\n")
    print(x$coefftable)

}

#' @export
summary.merlin <- function(object, ...) {

    cat("Mixed effects regression model \n")
    cat("Log likelihood = ")
    cat(object$loglikelihood)
    cat("\n")
    cat("\n")

    print(object$coefftable)

    if (length(object$levels)>0) {
        cat("\n")
        cat("Integration method: ")
        if (object$intmethod[1]=="ghermite")    cat("Non-adaptive Gauss-Hermite quadrature \n")
        else if (object$intmethod[1]=="halton") cat("Monte-Carlo integration using Halton sequences \n")
        else if (object$intmethod[1]=="sobol")  cat("Monte-Carlo integration using Sobol sequences \n")
        else                                    cat("Monte-Carlo integration \n")
        cat("Integration points: ")
        cat(object$ip)
        cat("\n")
        cat("       Convergence: ")
        cat(object$convergence)
    }

}

