library(simulator)
library(LEP)
library(IGESS)
library(GPA)
library(glmnet)
library(parallel)

#this block of codes could produce Figure
name_of_simulation <- "Along_U_Rho"
nrep <- 50
persize <- 2
sim <- new_simulation(name_of_simulation, "Compare Performance Along Different Pleiotropy Settings")
sim <-  generate_model(sim, make_sparse_linear_model,
                       n = 2000, p = 10000, h = 0.5, h0 = 0.5, s = 0.05, geno = 1, K=1, model_name = "slm_along_urho",
                       n0 = 8000,rho = as.list(c(0, 0.2, 0.4,0.6, 0.8)), u = as.list(c(0.9, 0.7, 0.3, 0.05)),vary_along = c("rho","u") )
nindex <- (nrep / persize)
for(index in 1:nindex){
  print(paste0("index:",index))
  sim <- simulate_from_model_parallel(sim, nsim = persize, index = index)
}
nthreads <- nindex
sim <-  run_method(sim, list(lep), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(igess), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(bvsr), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(gpa), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(Lasso), parallel = list(socket_names = nthreads, libraries = "glmnet"))
sim <- evaluate_parallel(sim, list(corr, rauc, power, sqrerr, FDR))
result <- as.data.frame(evals(sim))



#######################################################################
#this block of codes could produce Figure
name_of_simulation <- 'Along_K'
sim <- new_simulation(name_of_simulation, "Compare Performance Along K") %>%
  generate_model(LEPSimu::make_sparse_linear_model,
                 n = 2000, p = 10000, h = 0.5, model_name = "along_k",  h0 = 0.5, s = 0.05, geno = 1, u = 0.5,
                 n0 = 8000,rho = as.list(c(0,0.2,0.4,0.6,0.8)), K = as.list(c(1, 3, 5)),vary_along = c("rho","K") )
nrep <- 50
persize <- 2
nindex <- (nrep / persize)
for(index in 1:nindex){
  print(paste0("index:",index))
  sim <- simulate_from_model_parallel(sim, nsim = persize, index = index)
}
sim <-  run_method(sim, list(lep), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(igess), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(bvsr), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(gpa), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(Lasso), parallel = list(socket_names = nthreads, libraries = "glmnet"))
sim <- evaluate_parallel(sim, list(corr, rauc, power, sqrerr, FDR))
result <- as.data.frame(evals(sim))


#######################################################################
#this block of codes could produce Figure
name_of_simulation <- 'case_control_rho_K'
sim <- new_simulation(name_of_simulation, "case-control studies along K rho")
sim <- generate_model(sim, make_sparse_linear_model_cc,
                      n = 2000, p = 10000, h = 0.5, n0 = 8000, K=1, h0 = 0.5,s = 0.05, model_name = "cc_urho", geno = 1,
                      rho = as.list(c(0,0.2,0.4,0.6,0.8)), u = as.list(c(0.9, 0.7, 0.3, 0.05)), vary_along = c("rho","u"))
nrep <- 50
persize <- 2
nindex <- (nrep / persize)
nthreads <- nindex
for(index in 1:nindex){
  print(paste0("index:",index))
  sim <- simulate_from_model_parallel(sim, nsim = persize, index = index)
}
nthreads <- nindex
sim <-  run_method(sim, list(lep), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(igess), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(bvsr), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(gpa), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(LassoCC), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(Lasso), parallel = list(socket_names = nthreads, libraries = "glmnet"))
sim <- evaluate_parallel(sim, list(corr,cauc, rauc, power, sqrerr, FDR))
result <- as.data.frame(evals(sim))


#######################################################################
#this block of codes could produce Figure
name_of_simulation <- 'Assume_Dependent_pvalue'
sim <- new_simulation(name_of_simulation, "Assume Dependence for the pvalues")
sim <- generate_model(sim, make_sparse_linear_model_dependent_pvalue,
                      n = 2000, p = 10000, h = 0.5, n0 = 8000,model_name="dependent", u =0.9, v = 0.99, h0 = 0.5, s = 0.05, K = 2,
                      rho = as.list(c(0,0.2,0.4,0.6,0.8)), dr = as.list(seq(0.7,0.9,length = 3)), vary_along = c("rho","dr"))
nrep <- 50
persize <- 2
nindex <- (nrep / persize)
for(index in 1:nindex){
  print(paste0("index:",index))
  sim <- simulate_from_model_parallel(sim, nsim = persize, index = index)
}
nthreads <- nindex
sim <-  run_method(sim, list(lep), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(igess), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(bvsr), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(gpa), parallel = list(socket_names = nthreads))
sim <-  run_method(sim, list(Lasso), parallel = list(socket_names = nthreads, libraries = "glmnet"))
sim <- evaluate_parallel(sim, list(corr, rauc, power, sqrerr, FDR))
result <- as.data.frame(evals(sim))
