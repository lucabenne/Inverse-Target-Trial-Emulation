# Shared helpers for the revised GBSG2 ITTE pipeline.
# Baseline covariates stay on their original scale.
# Numerical standardisation is used only in the Weibull event/censoring models.

prepare_GBSG2_data <- function(){
  if(!requireNamespace("TH.data", quietly=TRUE)){
    stop("Package 'TH.data' is required.")
  }
  data("GBSG2", package="TH.data")
  datt <- GBSG2
  datt$tgrade <- as.numeric(datt$tgrade)
  datt$horTh <- as.integer(datt$horTh)-1L
  datt$menostat <- as.integer(datt$menostat)-1L
  colnames(datt)[1] <- "treat"
  colnames(datt)[10] <- "event"
  datt
}

GBSG2_survival_covariates <- c(
  "age","menostat","tsize","tgrade","pnodes","progrec","estrec"
)

get_GBSG2_survival_scaling <- function(data){
  means <- sapply(data[,GBSG2_survival_covariates,drop=FALSE], mean)
  sds <- sapply(data[,GBSG2_survival_covariates,drop=FALSE], sd)
  if(any(!is.finite(sds)) || any(sds<=0)){
    stop("All GBSG2 survival covariates must have finite positive SDs.")
  }
  list(mean=means,sd=sds)
}

apply_GBSG2_survival_scaling <- function(data, scaling){
  out <- data
  for(v in GBSG2_survival_covariates){
    out[[paste0(v,"_s")]] <- (out[[v]]-scaling$mean[[v]])/scaling$sd[[v]]
  }
  out
}

GBSG2_full_survival_formula <- function(){
  survival::Surv(time,event) ~ as.factor(treat) +
    age_s + menostat_s + tsize_s + tgrade_s + pnodes_s + progrec_s + estrec_s
}

# Linear predictor/scale for the fixed 9-coefficient full GBSG2 Weibull AFT DGP.
# survHEhmc parameterisation: T ~ Weibull(shape=alpha, scale=exp(X beta)).
GBSG2_weibull_scale <- function(data, beta, scaling, treat_override=NULL){
  if(length(beta)!=9L){
    stop("The full GBSG2 Weibull DGP requires exactly 9 beta coefficients.")
  }
  tr <- if(is.null(treat_override)) data$treat else rep(treat_override,nrow(data))
  z <- lapply(GBSG2_survival_covariates,function(v){
    (data[[v]]-scaling$mean[[v]])/scaling$sd[[v]]
  })
  names(z) <- GBSG2_survival_covariates
  eta <- beta[1] + beta[2]*tr +
    beta[3]*z$age + beta[4]*z$menostat + beta[5]*z$tsize +
    beta[6]*z$tgrade + beta[7]*z$pnodes + beta[8]*z$progrec + beta[9]*z$estrec
  sc <- exp(eta)
  if(any(!is.finite(sc)) || any(sc<=0)) stop("Invalid Weibull scale generated.")
  sc
}

GBSG2_draw_weibull_time <- function(data, beta, alpha, scaling, treat_override=NULL){
  if(!is.finite(alpha) || alpha<=0) stop("Weibull alpha must be positive.")
  sc <- GBSG2_weibull_scale(data,beta,scaling,treat_override)
  stats::rweibull(nrow(data),shape=alpha,scale=sc)
}
