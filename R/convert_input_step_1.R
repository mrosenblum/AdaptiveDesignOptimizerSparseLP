#' Adaptive Enrichment Design Optimization Using Sparse Linear Programming
#' Authors: Michael Rosenblum, Ethan Fang, Han Liu
#'
#' @param type.of.outcome.data "continuous", "binary"
#' @param subpopulation.1.size Proportion of overall population in subpopulation 1. Must be between 0 and 1.
#' @param total.alpha Familywise Type I error rate (1-sided)
#' @param accrual.yearly.rate Number of participants enrolled per year; assumed constant throughout trial
#' @param followup.length Time from enrollment to measurement of primary outcome (only used for continuous or binary outcome types)
#' @param population.parameters Matrix encoding scenarios (data generating distributions) used to define power constraints and  objective function; each row defines a data generating distribution; the columns represent subpopulation 1 outcome mean and variance followed by subpopulation 2 outcome mean and variance (4 entries total).
#' @param number.choices.end.of.stage.1 Number of enrollment choices at end of stage 1 
#' @param stage.1.sample.sizes Vector with 2 entries representing stage 1 sample sizes for subpopulations 1 and 2, respectively
#' @param stage.2.sample.sizes.per.enrollment.choice Matrix with number.choices.end.of.stage.1 rows and 2 columns, where the (i,j) entry represents the stage 2 sample size under enrollment choice i for subpopulation j.
#' @param objective.function.weights Vector with length equal to number of rows of population.parameters, representing weights used to define the objective function
#' @param power.constraints Matrix with same number of rows as population.parameters (each representing a data generating distribution) and three columns corresponding to the required power to reject (at least) H_01, H_02, H_0C, respectively.
#' @param type.of.LP.solver "matlab", "cplex", "GLPK" The linear program solve that you want to use; assumes that you have installed this already
#' @return 4 element list containing optimized designs from four classes (with increasing complexity):
#' @section Designs:
#' Each optimized design is a list containing: design.parameters and design.performance
#' @section design.parameters:
#' design.parameters has the following elements:
#' @export
#' @examples
#' #For demonstration purposes, the examples below use a coarse discretization.
#' #Example 1: continuous outcome
#' input_vector <- convert_input_step_1(
#'   ui.n.arms=2,
#'   ui.type.of.outcome.data="continuous",
#'   ui.time.to.event.trial.type="",
#'   ui.time.to.event.non.inferiority.trial.margin=NULL,
#'   ui.subpopulation.1.size=0.5,
#'   ui.total.alpha=0.05,
#'   ui.max.size=1000,
#'   ui.max.duration=5,
#'   ui.accrual.yearly.rate=250,
#'   ui.followup.length=1,
#'   ui.optimization.target="size",
#'   ui.time.to.event.censoring.rate=0,
#'   ui.mcid=NULL,
#'   ui.incorporate.precision.gain=FALSE,
#'   ui.relative.efficiency=1,
#'   ui.max.stages=5,
#'   ui.include.designs.start.subpop.1=FALSE,
#'   ui.population.parameters=matrix(c(15,15,3600,3600,3600,3600,
#'                                     15,0,3600,3600,3600,3600,
#'                                     0,15,3600,3600,3600,3600,
#'                                     0,0,3600,3600,3600,3600),nrow=4, ncol=6, byrow=TRUE,dimnames=
#'     list(c(),c("delta1","delta2","sigma1_trt","sigma1_con","sigma2_trt","sigma2_con"))),
#'   ui.desired.power=matrix(c(0,0,0.8,
#'                             0.8,0,0,
#'                             0,0.8,0,
#'                             0,0,0), nrow=4, ncol=3, byrow=TRUE,
#'     dimnames=list(c(),c("Pow_H(0,1)","Pow_H(0,2)","Pow_H(0,C)"))),
#'   ui.scenario.weights=matrix(c(0.25,0.25,0.25,0.25),ncol=1,dimnames=list(c(),c("weight"))),
#'   simulated.annealing.parameter.max.iterations=2
#' )
#'
#' #Example 2: binary outcome; 1 treatment arm versus control; superiority design
#' optimized_designs <- optimize_designs(
#'   ui.n.arms=2,
#'   ui.type.of.outcome.data="binary",
#'   ui.time.to.event.trial.type="",
#'   ui.time.to.event.non.inferiority.trial.margin=NULL,
#'   ui.subpopulation.1.size=0.4,
#'   ui.total.alpha=0.05,
#'   ui.max.size=2000,
#'   ui.max.duration=5,
#'   ui.accrual.yearly.rate=400,
#'   ui.followup.length=0,
#'   ui.optimization.target="size",
#'   ui.time.to.event.censoring.rate=0,
#'   ui.mcid=NULL,
#'   ui.incorporate.precision.gain=TRUE,
#'   ui.relative.efficiency=1.2,
#'   ui.max.stages=5,
#'   ui.include.designs.start.subpop.1=FALSE,
#'   ui.population.parameters=matrix(c(0.4,0.3,0.5,0.4,0.4,0.3,0.4,0.4,0.3,0.3,0.4,0.4),
#'     nrow=3, ncol=4, byrow=TRUE,dimnames=list(c(),c("p1_trt","p1_con","p2_trt","p2_con"))),
#'   ui.desired.power=matrix(c(0,0,0.8,0.8,0,0,0,0,0), nrow=3, ncol=3, byrow=TRUE,
#'     dimnames=list(c(),c("Pow_H(0,1)","Pow_H(0,2)","Pow_Reject_H0,1_and_H0,2"))),
#'   ui.scenario.weights=matrix(c(0.33,0.33,0.34),ncol=1,dimnames=list(c(),c("weight"))),
#'   simulated.annealing.parameter.max.iterations=2
#' )
#'
#'  #Example 3: continuous outcome; 2 treatment arms versus control; superiority design
#'  optimized_designs <- optimize_designs(
#'    ui.n.arms=3,
#'    ui.type.of.outcome.data="continuous",
#'    ui.time.to.event.trial.type="",
#'    ui.time.to.event.non.inferiority.trial.margin=NULL,
#'    ui.subpopulation.1.size=0.49,
#'    ui.total.alpha=0.05,
#'    ui.max.size=3000,
#'    ui.max.duration=8,
#'    ui.accrual.yearly.rate=240,
#'    ui.followup.length=0.5,
#'    ui.optimization.target="size",
#'    ui.time.to.event.censoring.rate=0,
#'    ui.mcid=NULL,
#'    ui.incorporate.precision.gain=TRUE,
#'    ui.relative.efficiency=1,
#'    ui.max.stages=4,
#'    ui.include.designs.start.subpop.1=FALSE,
#'    ui.population.parameters=matrix(c(0,3600,0,3600,0,3600,0,3600,0,3600,0,3600,0,3600,0,3600,
#'    15,3600,0,3600,0,3600,0,3600,0,3600,0,3600,15,3600,0,3600,15,3600,0,3600,0,3600,0,3600,15,
#'    3600,15,3600,0,3600,0,3600,0,3600,0,3600,15,3600,15,3600,15,3600,0,3600,0,3600,0,3600,15,
#'    3600,15,3600,15,3600,15,3600), nrow=6, ncol=12, byrow=TRUE),
#'    ui.desired.power=matrix(c(
#'      0,0,0,0,0,
#'      0.8,0,0,0,0,
#'      0.8,0,0.8,0,0,
#'      0.8,0.8,0,0,0,
#'      0.8,0.8,0.8,0,0,
#'      0.8,0.8,0.8,0.8,0),
#'      nrow=6, ncol=5, byrow=TRUE),
#'    ui.scenario.weights=matrix(c(0.166,0.166,0.166,0.166,0.166,0.167)),
#'    simulated.annealing.parameter.max.iterations=2)
#' @importFrom stats plogis
#' @importFrom mvtnorm pmvnorm GenzBretz
optimize_designs <- function(
  ui.n.arms,
  ui.type.of.outcome.data,
  ui.time.to.event.trial.type,
  ui.time.to.event.non.inferiority.trial.margin,
  ui.subpopulation.1.size,
  ui.total.alpha,
  ui.max.size,
  ui.max.duration,
  ui.accrual.yearly.rate,
  ui.followup.length,
  ui.optimization.target,
  ui.time.to.event.censoring.rate,
  ui.mcid,
  ui.incorporate.precision.gain,
  ui.relative.efficiency,
  ui.max.stages,
  ui.include.designs.start.subpop.1,
  ui.population.parameters,
  ui.desired.power,
  ui.scenario.weights,
  min.n.per.arm =25,       # For Continuous/Binary Outcomes
  min.enrollment.period =0.5,    # For Survival Outcomes
  simulated.annealing.parameter.function.scale =1,
  simulated.annealing.parameter.n.scale =100,
  simulated.annealing.parameter.period.scale =2,
  simulated.annealing.parameter.max.iterations =1000,
  simulated.annealing.parameter.n.simulations =1e4,
  simulated.annealing.parameter.means.temperature =100,
  simulated.annealing.parameter.survival.temperature =10,
  simulated.annealing.parameter.evals.per.temp =10,
  simulated.annealing.parameter.report.iteration =1,
  simulated.annealing.parameter.power.penalty =100000
){
  # Get start time
  isa.start.time <- proc.time()

  ### NOTE: restricted to two subpopulations ###
  n.subpopulations <- 2
  n.arms <- ui.n.arms
  ui.subpopulation.sizes <- c(ui.subpopulation.1.size, 1-ui.subpopulation.1.size)
  # If random seed is supplied, specify seeds. Otherwise pseudorandom seeds
  # are chosen based on the initial RNG state.
  if(!exists("initial.seed")){
    initial.seed <- sample(x=1:1e8, size=1)
  }

  # Set random seed
  set.seed(initial.seed)

  # Source design.evaluation code corresponding to number of arms in trial
  if(n.arms==2){
    # Computes distribution of test statistics in a given scenario,
    # using canonical joint distribution
    construct.joint.distribution.of.test.statistics <-
      function(...){
        construct.joint.distribution.of.test.statistics.OneTreatmentArm(...)
      }
    # Computes efficacy stopping boundaries
    generate.efficacy.boundaries <-
      function(...){
        get.eff.bound.OneTreatmentArm(...)
      }
    # Evaluates performance of simulated trials
    design.evaluate <-
      function(...){
        design.evaluate.OneTreatmentArm(...)
      }
    summarize.design.parameters.and.performance <-
      function(...){
        summarize.design.parameters.and.performance.OneTreatmentArm(...)
      }
  } 
  # Set functions for computing design features and design evaluation
  ##
  ## Format User Inputs from Graphical User Interface
  ##

  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    if(n.arms==2){
      if(ui.type.of.outcome.data=="binary") {
        ui.outcome.mean <- subset(ui.population.parameters,select=c(2,4,1,3))
        ui.outcome.sd <- sqrt(ui.outcome.mean*(1-ui.outcome.mean))
      } else{
        ui.outcome.mean <- cbind(array(0,c(nrow(ui.population.parameters),2)),subset(ui.population.parameters,select=c(1,2)))
        ui.outcome.sd <- sqrt(subset(ui.population.parameters,select=c(4,6,3,5)))
      }
    } else if(n.arms==3){
      if(ui.type.of.outcome.data=="binary") {
        ui.outcome.mean <- ui.population.parameters
        ui.outcome.sd <- sqrt(ui.outcome.mean*(1-ui.outcome.mean))
      } else{
        ui.outcome.mean <- subset(ui.population.parameters,select=c(1,3,5,7,9,11))
        ui.outcome.sd <- sqrt(subset(ui.population.parameters,select=c(2,4,6,8,10,12)))
      }
    }
    arm.names <- c(LETTERS[3], LETTERS[1:n.arms][-3])[1:n.arms]
    colnames(ui.outcome.sd) <- colnames(ui.outcome.mean) <-
      paste0(rep(arm.names, each=n.subpopulations),
             rep(1:n.subpopulations, n.arms))
  } 

  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    osoa.result <-
      sa.optimize(search.parameters=
                    list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual),
                         alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                              number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    # Cap sample size at minimum of the maximum specified size
                    # and the accrual rate x maximum allowable duration
                    list(n.per.arm=function(x)
                      ceiling(
                        squash(x,
                               min.n.per.arm,
                               min(ui.max.size, max.possible.accrual)/n.arms)),
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        delay=ui.followup.length,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type=ui.type.of.outcome.data,
                                        interim.info.times=NULL,
                                        outcome.mean=ui.outcome.mean,
                                        outcome.sd=ui.outcome.sd,
                                        mcid=ui.mcid,
                                        futility.boundaries=NULL,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.means.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)

  } else { # Survival Cases
    osoa.result <-
      sa.optimize(search.parameters=
                    list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                  min(c(feasible.max.duration,
                                                        ui.max.size/ui.accrual.yearly.rate))),
                         alpha.allocation=
                           rep(1/number.of.alpha.allocation.components,
                               number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    list(enrollment.period=function(x)
                      squash(x, min.enrollment.period,feasible.enrollment.period),
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type='survival',
                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                        hazard.rate=ui.hazard.rate,
                                        time=feasible.max.duration,
                                        max.follow=Inf,
                                        censoring.rate=ui.time.to.event.censoring.rate,
                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                        restrict.enrollment=FALSE,
                                        mcid=ui.mcid,
                                        futility.boundaries=NULL,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.period.scale,rep(1,number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.survival.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  }

  ## Two stage design
  n.stages <- 2 # Two Stage
  number.of.alpha.allocation.components <- n.stages*n.subpopulations

  ## 2SEA 2 stage equal alpha
  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    # if(ui.optimization.target=="ESS") {
    #   Switch Objective Function and Parameters
    # }
    tsea.result <-
      sa.optimize(search.parameters=
                    list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual)),
                  search.transforms=
                    # Cap sample size at minimum of the maximum specified size
                    # and the accrual rate x maximum allowable duration
                    list(n.per.arm=function(x)
                      ceiling(
                        squash(x,
                               min.n.per.arm,
                               min(ui.max.size, max.possible.accrual)/n.arms))
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        delay=ui.followup.length,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        interim.info.times=c(1/2,1),
                                        outcome.type=ui.type.of.outcome.data,
                                        outcome.mean=ui.outcome.mean,
                                        outcome.sd=ui.outcome.sd,
                                        mcid=ui.mcid,
                                        futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                             number.of.alpha.allocation.components),
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=simulated.annealing.parameter.n.scale,
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.means.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  } else { # Survival Cases
    if(ui.include.designs.start.subpop.1){
      number.of.alpha.allocation.components <- number.of.alpha.allocation.components - (n.subpopulations-1)}


    tsea.result <-
      sa.optimize(search.parameters=
                    list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                  min(feasible.max.duration,
                                                      ui.max.size/ui.accrual.yearly.rate))),
                  search.transforms=
                    list(enrollment.period=function(x)
                      squash(x, min.enrollment.period,
                             min(ui.max.duration,
                                 ui.max.size/ui.accrual.yearly.rate))),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type='survival',
                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                        hazard.rate=ui.hazard.rate,
                                        time=c(feasible.max.duration/2,feasible.max.duration),
                                        max.follow=Inf,
                                        censoring.rate=ui.time.to.event.censoring.rate,
                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                        restrict.enrollment=FALSE,
                                        mcid=ui.mcid,
                                        futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        alpha.allocation=
                                          rep(1/number.of.alpha.allocation.components,
                                              number.of.alpha.allocation.components),
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=simulated.annealing.parameter.period.scale,
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.survival.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  }

  ## 2SOA 2 stage optimized alpha
  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    tsoa.result <-
      sa.optimize(search.parameters=
                    list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual),
                         interim.info.times=c(1/2,1),
                         futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                         alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                              number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    # Cap sample size at minimum of the maximum specified size
                    # and the accrual rate x maximum allowable duration
                    list(n.per.arm=function(x)
                      ceiling(
                        squash(x,
                               min.n.per.arm,
                               min(ui.max.size, max.possible.accrual)/n.arms)),
                      interim.info.times=function(x){c(squash(x[1],0.1,0.9),1)},
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        delay=ui.followup.length,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type=ui.type.of.outcome.data,
                                        outcome.mean=ui.outcome.mean,
                                        outcome.sd=ui.outcome.sd,
                                        mcid=ui.mcid,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate
                  ),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.means.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)

  } else { # Survival Cases
    tsoa.result <-
      sa.optimize(search.parameters=
                    list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                  min(feasible.max.duration,
                                                      ui.max.size/ui.accrual.yearly.rate)),
                         time=c(feasible.max.duration/2,feasible.max.duration),
                         futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                         alpha.allocation=
                           rep(1/number.of.alpha.allocation.components,
                               number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    list(enrollment.period=function(x)
                      squash(x, min.enrollment.period,
                             min(ui.max.duration,
                                 ui.max.size/ui.accrual.yearly.rate)),
                      time=function(t){t1 <- squash(t[1],0.01,feasible.max.duration-0.02); t2<-squash(t[2],t1+0.01,feasible.max.duration); return(c(t1,t2))},
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type='survival',
                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                        hazard.rate=ui.hazard.rate,
                                        max.follow=Inf,
                                        censoring.rate=ui.time.to.event.censoring.rate,
                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                        restrict.enrollment=FALSE,
                                        mcid=ui.mcid,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.survival.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  }

  ####
  # Optimize 2 Stage, Group Sequential Design, for 2 arm trials
  ####
  if(n.arms==2){
    # Computes distribution of test statistics in a given scenario,
    # using canonical joint distribution
    construct.joint.distribution.of.test.statistics <-
      function(...){
        construct.joint.distribution.of.test.statistics.GroupSequential.OneTreatmentArm(...)
      }
    # Computes efficacy stopping boundaries
    generate.efficacy.boundaries <-
      function(...){
        get.eff.bound.GroupSequential.OneTreatmentArm(...)
      }
    # Evaluates performance of simulated trials
    design.evaluate <-
      function(...){
        design.evaluate.GroupSequential.OneTreatmentArm(...)
      }
    ## Optimize:
    if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
      group.sequential.tsoa.result <-
        sa.optimize(search.parameters=
                      list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual),
                           interim.info.times=c(1/2,1),
                           futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                           alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                number.of.alpha.allocation.components)
                      ),
                    search.transforms=
                      # Cap sample size at minimum of the maximum specified size
                      # and the accrual rate x maximum allowable duration
                      list(n.per.arm=function(x)
                        ceiling(
                          squash(x,
                                 min.n.per.arm,
                                 min(ui.max.size, max.possible.accrual)/n.arms)),
                        interim.info.times=function(x){c(squash(x[1],0.1,0.9),1)},
                        alpha.allocation=reals.to.probability
                      ),
                    fixed.parameters=list(n.arms=n.arms,
                                          accrual.rate=ui.accrual.yearly.rate,
                                          delay=ui.followup.length,
                                          subpopulation.sizes=ui.subpopulation.sizes,
                                          outcome.type=ui.type.of.outcome.data,
                                          outcome.mean=ui.outcome.mean,
                                          outcome.sd=ui.outcome.sd,
                                          mcid=ui.mcid,
                                          relative.efficiency=ui.relative.efficiency,
                                          n.simulations=simulated.annealing.parameter.n.simulations,
                                          total.alpha=ui.total.alpha,
                                          construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                          generate.efficacy.boundaries=generate.efficacy.boundaries,
                                          design.evaluate=design.evaluate
                    ),
                    create.object=triage.based.on.outcome.type,
                    evaluate.object=power.penalized.weighted,
                    function.scale=simulated.annealing.parameter.function.scale,
                    parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                    max.iterations=simulated.annealing.parameter.max.iterations,
                    temperature=simulated.annealing.parameter.means.temperature,
                    evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                    report.iteration=simulated.annealing.parameter.report.iteration,
                    scenario.weights=ui.scenario.weights,
                    power.penalty=simulated.annealing.parameter.power.penalty,
                    power.constraints=ui.desired.power,
                    optimization.target=ui.optimization.target)

    } else { # Survival Cases
      group.sequential.tsoa.result <-
        sa.optimize(search.parameters=
                      list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                    min(feasible.max.duration,
                                                        ui.max.size/ui.accrual.yearly.rate)),
                           time=c(feasible.max.duration/2,feasible.max.duration),
                           futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                           alpha.allocation=
                             rep(1/number.of.alpha.allocation.components,
                                 number.of.alpha.allocation.components)
                      ),
                    search.transforms=
                      list(enrollment.period=function(x)
                        squash(x, min.enrollment.period,
                               min(ui.max.duration,
                                   ui.max.size/ui.accrual.yearly.rate)),
                        time=function(t){t1 <- squash(t[1],0.01,feasible.max.duration-0.02); t2<-squash(t[2],t1+0.01,feasible.max.duration); return(c(t1,t2))},
                        alpha.allocation=reals.to.probability
                      ),
                    fixed.parameters=list(n.arms=n.arms,
                                          accrual.rate=ui.accrual.yearly.rate,
                                          subpopulation.sizes=ui.subpopulation.sizes,
                                          outcome.type='survival',
                                          non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                          hazard.rate=ui.hazard.rate,
                                          max.follow=Inf,
                                          censoring.rate=ui.time.to.event.censoring.rate,
                                          ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                          restrict.enrollment=FALSE,
                                          mcid=ui.mcid,
                                          relative.efficiency=ui.relative.efficiency,
                                          n.simulations=simulated.annealing.parameter.n.simulations,
                                          total.alpha=ui.total.alpha,
                                          construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                          generate.efficacy.boundaries=generate.efficacy.boundaries,
                                          design.evaluate=design.evaluate),
                    create.object=triage.based.on.outcome.type,
                    evaluate.object=power.penalized.weighted,
                    function.scale=simulated.annealing.parameter.function.scale,
                    parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                    max.iterations=simulated.annealing.parameter.max.iterations,
                    temperature=simulated.annealing.parameter.survival.temperature,
                    evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                    report.iteration=simulated.annealing.parameter.report.iteration,
                    scenario.weights=ui.scenario.weights,
                    power.penalty=simulated.annealing.parameter.power.penalty,
                    power.constraints=ui.desired.power,
                    optimization.target=ui.optimization.target)
    }
  }
  if(ui.n.arms==2){
    output_summary <- lapply(list(Single.Stage.Equal.Alpha.Allocation.Design=osea.result,Single.Stage.Optimized.Alpha.Allocation.Design=osoa.result,Two.Stage.Group.Sequential.Design=group.sequential.tsoa.result,Two.Stage.Equal.Alpha.Allocation.Design=tsea.result,Two.Stage.Optimized.Alpha.Allocation.Design=tsoa.result),summarize.design.parameters.and.performance)}else{
      output_summary <- lapply(list(Single.Stage.Equal.Alpha.Allocation.Design=osea.result,Single.Stage.Optimized.Alpha.Allocation.Design=osoa.result,Two.Stage.Equal.Alpha.Allocation.Design=tsea.result,Two.Stage.Optimized.Alpha.Allocation.Design=tsoa.result),summarize.design.parameters.and.performance)}
  return(output_summary)
}
