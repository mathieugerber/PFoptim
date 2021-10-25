
#' Global Particle filter Stochastic Optimization
#'
#' @importFrom mnormt  rmt  rmnorm 
#' @importFrom Rdpack reprompt
#' @description This function implements the G-PFSO (Global Particle Filter Stochastic Optimization) algorithm of \insertCite{gerber2020online2;textual}{PFoptim} for minimzing either the function \eqn{\theta\mapsto E[\mathrm{fn}(\theta,Y)]} from i.i.d. realizations \eqn{y_1,...,y_n} of \eqn{Y} or the function  \eqn{\theta\mapsto\sum_{i=1}^n \mathrm{fn}(\theta,y_i)}, where \eqn{\theta} is a vector of dimension d.

#' @usage gpfso(y, N, fn, init, numit, ..., resampling=c("SSP", "STRAT", "MULTI"), control= list())
#' @param y Either a vector of observations or a matrix of observations (the number of rows being the sample size). 
#' @param N Number of particles. The parameter  \code{N}  must be greater or equal to 2.
#' @param fn function for a single observation. If theta is an  \code{N} by d matrix and \code{y} is a matrix then  \code{fn(theta, y[i,])} must be a vector of length \code{N}. Similarly, if theta is an  \code{N} by d matrix  and \code{y} is a vector then \code{fn(theta,y[i])} must be a vector of length \code{N}. If some rows of theta are outside the search space then the corresponding entries of the vector \code{fn(theta, y[i,])} must be equal to \code{Inf}.
#' @param init Function used to sample the initial particles  such that \code{init(N)} is an \code{N} by d matrix (or alternatively a vector of length \code{N} if d=1).
#' @param ... Further arguments to be passed to \code{fn}.
#' @param numit Number of iterations of the algorithm. If \code{numit} is not specified then G-PFSO estimates the minimizer of the function \eqn{E[\mathrm{fn}(\theta,Y)]} (in which case the observations are processed sequentially and \code{numit} is equal to the sample size). If \code{numit} is specified then G-PFSO computes the minimizer of the function \eqn{\sum_{i=1}^n \mathrm{fn}(\theta,y_i)}.
#' @param resampling Resampling algorithm to be used. Resamping should be either "SSP" (SSP resampling), "STRAT" (stratified resampling) or "MULTI" (multinomial resampling).
#' @param control A \code{list} of control parameters. See details.
#' @details 
#' Note that arguments after \code{...} must be matched exactly.
#'
#' G-PFSO computes two estimators of the minimizer of the objective function, namely the estimators \eqn{\bar{\theta}^N_{\mathrm{numit}}} and  \eqn{\tilde{\theta}^N_{\mathrm{numit}}}. The former  is defined by \eqn{\bar{\theta}^N_{\mathrm{numit}}=\frac{1}{\mathrm{numit}}\sum_{t=1}^{\mathrm{numit}}\tilde{\theta}_t^N}  and converges to a particular element of the search space at a faster rate than the latter, but the latter estimator can find more quickly a small neighborhood of the minimizer of the objective function. 
#'
#' By default the sequence \eqn{(t_p)_{p\ge  0}} is taken as 
#'\deqn{t_p=t_{p-1}+\lceil  \max\big( A t_{p-1}^{\varrho}\log(t_{p-1}),B\big) \rceil}
#' with A=B=1, \eqn{\varrho=0.1} and \eqn{t_0=5}. The value of \eqn{A,B,\varrho} and \eqn{t_0} can be changed using the control argument (see below).
#'
#' The \code{control} argument is a list that can supply any of the following components:
#' \describe{
#' \item{alpha:}{Parameter \eqn{\alpha} of the learning rate \eqn{t^{-\alpha}}, which must be a strictly positive real number. By default,  \code{alpha=0.5}.}
#' \item{Sigma:}{Scale matrix used to sample the particles.  \code{Sigma} must be either a d by d covariance matrix or a strictly positive real number. In this latter case the scale matrix used to sample the particles is  \code{diag(Sigma , d )}. By default,  \code{Sigma=1}.}
#' \item{trace:}{If trace=TRUE then the value of \eqn{\tilde{\theta}_{t}} and of the effective sample size \eqn{ESS_t} for all \eqn{t=1,\dots,\mathrm{numit}} are returned. By default, trace=FALSE.}
#' \item{A:}{Parameter A of the sequence \eqn{(t_p)_{p\geq 0}} used by default (see above). This parameter must be strictly positive.}
#' \item{B:}{Parameter B of the sequence \eqn{(t_p)_{p\geq 0}} used by default (see above). This parameter must non-negative.}
#' \item{varrho:}{Parameter varrho of the sequence \eqn{(t_p)_{p\geq 0}} used by default (see above). This parameter must be in the interval (0,1).}
#' \item{t0:}{Parameter \eqn{t_0} of the the sequence \eqn{(t_p)_{p\geq 0}} used by default (see above). This parameter must be a non-negative integer.}
#' \item{nu:}{Number of degrees of freedom of the Student's t-distributions used at time \eqn{t\in(t_p)_{\ge 0}} to generate the new particles. By default \code{nu=10}}
#' \item{c_ess:}{A resamling step is performed when \eqn{ESS_t<=N c_\mathrm{ess}}. This parameter must be in the interval (0,1] and by default \code{c_ess}=0.7.}
#'  }
#' @return A list with the following components:
#'   \item{B_par}{Value of \eqn{\bar{\theta}^N_{\mathrm{numit}}}}
#'   \item{T_par}{Value of \eqn{\tilde{\theta}^N_{\mathrm{numit}}}}
#'   \item{T_hist}{Value of \eqn{\tilde{\theta}^N_{t}} for \eqn{t=1,...,\mathrm{numit}} (only if trace=TRUE)}
#'   \item{ESS}{Value of the effective sample for \eqn{t=1,...,\mathrm{numit}} (only if trace=TRUE)}
#' @keywords global optimization, stochastic optimization, particle filters
#' @export
#' @examples
#' #Definition of fn
#' fn_toy<-function(theta, obs){
#'   test<-rep(0,nrow(theta))
#'   test[theta[,2]>0]<-1
#'   ll<-rep(-Inf,nrow(theta))
#'   ll[test==1]<-dnorm(obs,mean=theta[test==1,1], sd=theta[test==1,2],log=TRUE)
#'   return(-ll)
#' }
#' #Generate data y_1,...,y_n
#' n<-10000             #sample size
#' theta_star<-c(0,1)  #true parameter value
#' y<-rnorm(n,mean=theta_star[1], sd=theta_star[2])
#' d<-length(theta_star)
#' #Define init funciton to be used
#' pi0<-function(N){
#'     return(cbind(rnorm(N,0,5), rexp(N)))
#' }
#' ##Example 1: Maximum likelihood estimation in the Gaussian model 
#' ##true value of the MLE
#' mle<-c(mean(y),sd(y))
#' ## use gpfso to compute the MLE
#' Est<-gpfso(y, N=100, fn=fn_toy, init=pi0, numit=20000, control=list(trace=TRUE))
#' ## print \bar{\theta}^N_{numit} and \tilde{\theta}^N_{numit}
#' print(Est$B_par)
#' print(Est$T_par)
#' ##assess convergence
#' par(mfrow=c(1,2))
#' for(k in 1:2){
#'   plot(Est$T_hist[,k],type='l', xlab="iteration", ylab="approximation error")
#'   lines(cumsum(Est$T_hist[,k])/1:length(Est$T_hist[,k]),type='l', col='red')
#'   abline(h=mle[k])
#' }
#' ##Example 2: Expected log-likelihood estimation in the Gaussian model 
#' ## Estimation of theta_star using gpfso
#' Est<-gpfso(y, N=100, fn=fn_toy, init=pi0, control=list(trace=TRUE))
#' ## print \bar{\theta}^N_{numit} and \tilde{\theta}^N_{numit}
#' print(Est$B_par)
#' print(Est$T_par)
#' ##assess convergence
#' par(mfrow=c(1,2))
#' for(k in 1:2){
#'   plot(Est$T_hist[,k],type='l', xlab="iteration", ylab="approximation error")
#'   lines(cumsum(Est$T_hist[,k])/1:length(Est$T_hist[,k]),type='l', col='red')
#'   abline(h=theta_star[k])
#' }
#' @references 
#' \insertAllCited{}

gpfso<-function(y, N, fn, init, ..., numit=-1, resampling="SSP", control= list()){
  #default values
  alpha<-0.5
  A<-1
  B<-1
  varrho<-0.1
  t0<-5
  c_ess<-0.7
  nu<-10
  
  if(is.numeric(N)==FALSE || length(c(N))!=1 || N<2 || N%%1!=0){
        return(cat('Error: N should be an integer greater or equal to 2'))
  }
  if(is.function(init)==FALSE || is.null(dim(init(N)))==TRUE){
       return(cat('Error: init should be a function such that init(N) is an N by d matrix (d= dimention of the search space) or a vector of length N (if d=1)'))
  }
  test<-init(2)
  if(is.vector(test)==TRUE){
        d<-1
        Sigma_use<-diag(1,d)
  }else if(is.matrix(test)==TRUE){
        d<-ncol(test)
        Sigma_use<-diag(1,d)
  }else{
        return(cat('Error: init should be a function such that init(N) is an N by d matrix (d= dimention of the search space) or a vector of length N (if d=1)'))
  }  
  if(is.numeric(numit)==FALSE || length(c(numit))!=1 || (numit<2 && numit!= -1) || numit%%1!=0){
        return(cat('Error: numit should be an integer greater or equal to 2'))
  }
  if(is.null(control$alpha)==FALSE){
        if(is.numeric(control$alpha)==FALSE || length(c(control$alpha))!=1 || control$alpha<=0){
              return(cat('Error: alpha should be a strictly positive real number'))
        }else{
              alpha<-control$alpha
        }
  }
  
  if(is.null(control$Sigma)==FALSE){
        if(is.numeric(control$Sigma)==FALSE){
               return(cat('Error: Sigma should be a either a stictly positive real number or a symmetric positive definite matrix')) 
        }else if(length(c(control$Sigma))==1 && control$Sigma<=0){
               return(cat('Error: Sigma should be a either a stictly positive real number or a symmetric positive definite matrix'))   
        }else if(length(c(control$Sigma))!=1 && is.matrix(control$Sigma)==FALSE){
               return(cat('Error: Sigma should be a either a stictly positive real number or a symmetric positive definite matrix'))        
        }else if( (length(c(control$Sigma))!=1 && is.matrix(control$Sigma)) && 
               (nrow(control$Sigma)!=d || ncol(control$Sigma)!=d   || max(abs(control$Sigma-t(control$Sigma)))>0 || min(eigen(control$Sigma)$values)<=0) ) {
                return(cat('Error: Sigma should be a either a stictly positive real number or a symmetric positive definite matrix'))  
        }
        if(length(c(control$Sigma))==1){
            Sigma_use<-diag(control$Sigma,d)
            chol_Sig<-chol(Sigma_use)
        }else{
            Sigma_use<-control$Sigma
            chol_Sig<-chol(Sigma_use)
        }
  }
  
  if(is.null(control$A)==FALSE){
           if(is.numeric(control$A)==FALSE || length(c(control$A))>1 || control$A<=0){
               return(cat('Error: A should be strictly positive real number'))
           }else{
               A<-control$A  
           }
  }
  if(is.null(control$B)==FALSE){
           if(is.numeric(control$B)==FALSE || length(c(control$B))>1 || control$B<0){
               return(cat('Error: B should be a non-negative real number'))
           }else{
               B<-control$B  
           }
  }
  if(is.null(control$varrho)==FALSE){
           if(is.numeric(control$varrho)==FALSE || length(c(control$varrho))>1 || control$varrho<=0 ||  control$varrho>=1){
              return(cat('Error: varrho should be a real number in (0,1)'))
           }else{
              varrho<-control$varrho
              if(varrho>=alpha){
                 cat('Warning: varrho should be smaller than alpha (by default alpha=0.5)')
              }
           }
  }        
  if(is.null(control$c_ess)==FALSE){
           if(is.numeric(control$c_ess)==FALSE || length(c(control$c_ess))>1 || control$c_ess<=0 || control$c_ess>1){
              return(cat('Error: c_ess should be in (0,1]'))
           }else{
              c_ess<-control$c_ess
           }
  }   
  if(is.null(control$t0)==FALSE){
           if(is.numeric(control$t0)==FALSE || length(c(control$t0))>1 || control$t0<0 || control$t0%%1!=0){
              return(cat('Error: t0 should be non-negative integer'))
           }else{
              t0<-control$t0
           } 
  }      
  if(is.null(control$nu)==FALSE){
           if(is.numeric(control$nu)==FALSE || length(c(control$nu))>1 || control$nu<=1){
              return(cat('Error: nu should be a scalar strictly greater than 1'))
           }else{
              nu<-control$nu
           }
    }
  if(is.numeric(y)==FALSE){
        return(cat('Error: y should be a numeric vector or a numeric matrix'))
  }else if(is.vector(y)==FALSE && is.matrix(y)==FALSE){
        return(cat('Error: y should be a numeric vector or a numeric matrix'))   
  }else if((resampling %in% c("SSP", "STRAT", "MULTI"))==FALSE){
       return(cat('Error: resampling should be either "SSP", or "STRAT", or "MULTI" '))
  }else{
       if(is.vector(y)){
          nobs<-length(y)
       }else{
          nobs<-nrow(y)
       }
       if(numit==-1){
         nit<-nobs
       }else{
         nit<- numit
       }
       
       if(is.null(control$tp)){  
         tp_use<-default_Tp(alpha,nit, A, B, varrho, t0)
       }else{
         tp_use<-c(control$tp)
       }
       
       if(resampling=="SSP"){
          Algo_res<-SSP_Resampler
       }else if(resampling=="STRAT"){
          Algo_res<-Stratified_Resampler
       }else{
          Algo_res<-function(U,W) return(sample(1:length(W),length(W),prob=W, replace=TRUE))
       }
       
       collapse<-FALSE
       ESS_bound<- N*c_ess 
       w<-rep(0,N)
       ESS_vec<-0  
       mean_vec<-0 
       if(is.null(control$trace)==FALSE && control$trace==TRUE){
          mean_vec<- matrix(0,nit,d)
          ESS_vec<-  rep(0,nit)
       }
       #initialization
       t<-1
       particles<-as.matrix(init(N))  
       if(numit==-1){
          use<-t
       }else{
          use<-sample(1:nobs,1)
       }
       if(is.vector(y)){
           work<- -fn(particles,   y[use], ...)
       }else{
           work<- -fn(particles,   y[use,],...)
       }
       if(is.numeric(work)==FALSE || length(c(work))!=N){
              return(cat('Error: for a single observation y, fn(particles, y) should be a vector of length N'))
       }else{
         work[is.nan(work)]<--Inf
         if(max(work)== -Inf){
              return(cat('Error: none of the particles returned by init(N) is in the search space'))
         }else{
           w<-work
           w1<- exp(w - max(w))
           W<- w1 / sum(w1)
           ESS<-1/sum(W^2)
           tilde_theta<-apply(W*particles,2,sum)
           bar_theta<-tilde_theta
           if(is.null(control$trace)==FALSE && control$trace==TRUE){
             mean_vec[t,]<-tilde_theta
             ESS_vec[t]<-ESS
           }
           count=1
           for(t in  2:nit){
            A<-1:N
            if(ESS<=ESS_bound){
               A<-Algo_res(runif(N),W)
               w<-rep(0,N)
            }
            if(t-1==tp_use[count]){
               particles<-as.matrix(particles[A,])+rmt(N, rep(0,d), Sigma_use*(t-1)^(-2*alpha), df=nu, sqrt=chol_Sig*(t-1)^(-alpha))
               count<-count+1
            }else{
               particles<-as.matrix(particles[A,])+rmnorm(N,rep(0,d), Sigma_use*(t-1)^(-2*alpha), sqrt=chol_Sig*(t-1)^(-alpha))
            } 
            if(numit==-1){
                use<-t
            }else{
                use<-sample(1:nobs,1)
            }
            if(is.vector(y)){
               work<- -fn(particles,   y[use],...)
            }else{
               work<- -fn(particles,   y[use,],...)
            }
            work[is.nan(work)]<--Inf
            if(max(work)==-Inf){
                collapse<- TRUE
                break
            }
            w<-w+work
            w1<- exp(w - max(w))
            W<- w1 / sum(w1)
            ESS<-1/sum(W^2)
            tilde_theta<-apply(W*particles,2,sum)
            bar_theta<-((t-1)*bar_theta+tilde_theta)/t
            if(is.null(control$trace)==FALSE && control$trace==TRUE){
               mean_vec[t,]<-tilde_theta
               ESS_vec[t]<-ESS
            }
        }
       }
     }
     if(collapse){
          return(cat('Error: particle system collapse (all the particles are outside the search space) '))
     }else{
           if(is.null(control$trace)==FALSE && control$trace==TRUE){
              return(list(B_par=bar_theta, T_par=tilde_theta, T_hist=mean_vec,  ESS=ESS_vec))
           }else{
              return(list(B_par=bar_theta, T_par=tilde_theta, T_hist=NULL,  ESS=NULL))
           }
     }
  }
}


 












