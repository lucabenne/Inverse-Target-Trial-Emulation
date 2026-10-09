library(R2jags)

generate_BG = function(data){

filein = '../src/04_GBSG2_covariate_DGP_model.txt'
  
  datalist=list(
    N = nrow(data),
    age = data$age,
    menostat = data$menostat,
    tsize = data$tsize,
    tgrade = as.integer(data$tgrade),
    pnodes = data$pnodes,
    sqrt_progrec = sqrt(data$progrec),
    sqrt_estrec = sqrt(data$estrec)
  )

  params = c(
    'theta0', 'theta10', 'theta11',
    'theta2', 'theta21', 'theta22',
    'grade_cut1', 'grade_delta', 'grade_cut2',
    'theta31', 'theta32', 'theta33',
    'theta40', 'theta41', 'theta42', 'theta43', 'theta44',
    'theta5', 'theta51', 'theta52', 'theta53', 'theta54', 'theta55',
    'theta6', 'theta61', 'theta62', 'theta63', 'theta64', 'theta65', 'theta66',
    'sigma0', 'sigma2', 'sigma4', 'sigma5', 'sigma6'
  )

  n.iter = 1000
  inits = NULL

  BG = jags(
    data=datalist,
    inits=inits,
    parameters.to.save=params,
    model.file=filein,
    n.chains=2,
    n.burnin=0,
    n.iter=n.iter,
    n.thin=1,
    DIC=FALSE
  )

  return(BG)
}
