default_Tp<-function(alpha, Nit, A, B, varrho, t0){
   Tp<-t0
   t1<-0
   while(t1 <= Nit){
     t1<-t0+ceiling(max(A*t0^varrho*log(t0),B))
     if(t1 <= Nit) Tp<-c(Tp, t1)
     t0<-t1
  }
  return(c(Tp,Nit+1))
}
