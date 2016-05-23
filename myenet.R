source("enet_funcs.R")

myarrayspc<-
function(x,K=1,para,use.corr=FALSE, max.iter=100,trace=FALSE,eps=1e-3)
{     call <- match.call()
      x<-scale(x,center=TRUE,scale=use.corr)
      svdobj<-svd(x)
      v<-svdobj$v
      totalvariance<-sum((svdobj$d)^2)      
      alpha<-as.matrix(v[,1:K,drop=FALSE])      
      beta<-alpha      
      for ( i in 1:K) {
         y<-drop(x%*%alpha[,i])
         beta[,i]<-soft(drop(t(x)%*%y),para=para[i])
         }
      temp<-beta
      normtemp<-sqrt(apply(temp^2,2,sum))
      normtemp[normtemp==0]<-1
      temp<-t(t(temp)/normtemp)
      k<-0
      diff<-1
      while((k<max.iter) & (diff>eps)){
         k<-k+1
         alpha<-x%*%beta
         alpha<-t(x)%*%alpha
         z<-svd(alpha)
         alpha<-(z$u)%*%t(z$v)
         for ( i in 1:K) {
           y<-drop(x%*%alpha[,i])
           beta[,i]<-soft(drop(t(x)%*%y),para=para[i])
         }     
         normbeta<-sqrt(apply(beta^2,2,sum))
         normbeta[normbeta==0]<-1
         beta2<-t(t(beta)/normbeta)
         diff<-max(abs(beta2-temp))
         temp<-beta2
          if(trace){
              if (k%%10==0){
                  cat("iterations",k,fill=TRUE)
                 }
            }
      } 
      normbeta<-sqrt(apply(beta^2,2,sum))
      normbeta[normbeta==0]<-1
      beta<-t(t(beta)/normbeta)
      u<-x%*%beta
      R<-qr.R(qr(u))
      pev<-diag(R^2)/totalvariance
      obj<-list(call = call, K=K,loadings=beta,pev=pev,var.all=totalvariance,para=para)
      class(obj) <- "arrayspc"
      obj 
}


