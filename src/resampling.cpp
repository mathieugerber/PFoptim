
#include <Rcpp.h>
using namespace Rcpp;


//' SSP resampling
//'
//' @importFrom Rdpack reprompt
//' @description This function implements the SSP resampling algorithm \insertCite{gerber2019negative}{PFoptim}.
//' @usage SSP_Resampler(U,W)
//' @param U A vector of points in (0,1) such that \code{length(U)=length(W)}.
//' @param W A vector of normalized weights.
//' @details 
//' For efficiency reasons, \code{SSP_Resampler} does not perform checks on the supplied arguments.
//' @return A vector of length N with elements in the set \eqn{\{1,...,N\}}, with \code{N=length(U)=length(W)}.
//' @keywords resampling algorithms
//' @export
//' @examples
//' N<-100
//' W<-rbeta(N,0.5,2)
//' W<-W/sum(W)
//' J<-SSP_Resampler(runif(N),W)
//' @references 
//' \insertAllCited{}
// [[Rcpp::export(SSP_Resampler)]]
 NumericVector SSP_Resampler(const NumericVector points,  const NumericVector W){
 
	int i,j,k, index1, index2;
	double s,p1,p2;
	int N=points.size();
	NumericVector J(N);

	j=0;
	for(i=0;i<N;i++)
	{
	       for(k=0;k<floor(N*W[i]); k++)
	       {
               	J[j]=i+1;
			j++;
		}
	}

        if(j<N){
	   s=0;
	   index1=0;
	   index2=1;
	   p1=N*W[index1]-floor(N*W[index1]);
	   for(i=0;i< N;i++)
	   {
	        p2=N*W[index2]-floor(N*W[index2]);
		s=p1+p2;
		if(s<1.0)
		{
			if(points[i]<p1/s)
			{
				p1=s;
				index2++;
			}
			else
			{
			        index1=index2;
				 p1=s;
				 index2++;
			}
		}
		else if(s==1)
		{
			if(points[i]<p1)
			{
				J[j]=index1+1;
				index1=index2+1;
				index2=index2+2;
				p1=N*W[index1]-floor(N*W[index1]);
				j++;
			}
			else
			{
				J[j]=index2+1;
				index1=index2+1;
				index2=index2+2;
				p1=N*W[index1]-floor(N*W[index1]);
				j++;
			}
		}
		else
		{
			if(points[i]<(1.0-p1)/(2.0-s))
			{
				p1=s-1.0;
				J[j]=index2+1;
				index2++;
				j++;
			}
			else
			{
			        J[j]=index1+1;
			        index1=index2;
				p1=s-1.0;
				index2++;
				j++;
			}
		}
		if(index2==N)
		{
			break;
		}
	   }
	   if(j==N-1)
	   {
		  J[N-1]=index1+1;
	   }
	}
	return J;
}



//' Stratified resampling
//'
//' @importFrom Rdpack reprompt
//' @description This function implements the stratified resampling algorithm descibed see e.g. in Section 9.6 of \insertCite{chopin2020introduction;textual}{PFoptim}
//' @usage Stratified_Resampler(U,W)
//' @param U A vector of points in (0,1) such that \code{length(U)=length(W)}.
//' @param W A vector of normalized weights.
//' @details 
//' For efficiency reasons, \code{Stratified_Resampler} does not perform checks on the supplied arguments.
//' @return A vector of length N with elements in the set \eqn{\{1,...,N\}}, with \code{N=length(U)=length(W)}.
//' @keywords resampling algorithms
//' @export
//' @examples
//' N<-100
//' W<-rbeta(N,0.5,2)
//' W<-W/sum(W)
//' J<-Stratified_Resampler(runif(N),W)
//' @references 
//' \insertAllCited{}
// [[Rcpp::export(Stratified_Resampler)]]
 NumericVector Stratified_Resampler(const NumericVector points,  const NumericVector W){
 
	int i,j,N=points.size();
	double s;
	NumericVector J(N);
	j=0;
	s=0;
	for(i=0;i<N;i++)
	{
		s+= W[i];
		while((j+points[j])/N<=s && j<N)
		{
			J[j]=i+1;
			j++;
		}

		if(j==N)
		{
			break;
		}
	}
	return(J);
}


















