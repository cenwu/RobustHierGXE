Data1 <- function(n,s)
{

         q = 10
         sige = matrix(0,q,q)

         for (i in 1: q)
            {
                   for (j in 1: q)
                    {
                      sige[i,j] = 0.5^abs(i-j)     # AR structure 
                    }           

            }
    e = mvrnorm(n,rep(0,q),sige)
    e = scale(e)  
    sig=matrix(0,s,s)

    # simulate the genetic factor 
         for (i in 1: s)
            {
                   for (j in 1: s)
                    {
                        sig[i,j] = 0.5^abs(i-j)                                        # AR  structure
                       # sig[i,j] = 1*(i==j) + 0.5*(abs(i-j)==1) + 0.33*(abs(i-j)==2)  # Banded  structure
                    }           

            }

    x = mvrnorm(n,rep(0,s),sig)

 #   SNP data
 #   for(j in 1:s)
 #          {
 #             t = x[,j]
 #             d = quantile(t,probs=c(1/3,2/3))
 #             names(d) = NULL
 #             x[,j] = I(t>d[2])*1 + I(t>d[1])*1
 #          }

    x = scale(x)  

    # x = cbind(1,x)
    #  err = rnorm(n,0,1)   
    # err = as.matrix(c(rnorm(0.80*n), 1*rlnorm(0.20*n)))
      err = 1*rt(n,df=2)
    # err = as.matrix(c(rnorm(0.80*n), 1*rcauchy(0.20*n)))
     err = err[sample(n)]

####################################################################################
 
     coef = runif(35,0.6,0.9)
     coef[1:q] = runif(10,0.8,0.8)
     coef[(q+1):(q+10)] = runif(10,0.8,1.2) 
     coef[(q+11):(q+25)] = runif(15,0.8,1.2) 

      mat = cbind(e,x[,(1:10)],scale(cbind(x[,1]*e[,1],x[,1]*e[,3],x[,1]*e[,5],x[,2]*e[,1],x[,2]*e[,3],x[,2]*e[,5],
                                           x[,3]*e[,1],x[,3]*e[,3],x[,3]*e[,5],x[,4]*e[,1],x[,4]*e[,3],x[,4]*e[,5],
                                           x[,5]*e[,1],x[,5]*e[,3],x[,5]*e[,5] )  ))


     y = rowSums(sweep(mat,MARGIN=2,coef,`*`)) + err   
     y1 = y
     LOG = min(c(y1))        
      tt = exp(c(y1))        
     tt[which(tt=="NaN"|tt=="Inf")] =  max(tt[which(tt!="NaN"&tt!="Inf")])
     m1 = 0.5*quantile(tt, 0.77)
     m2 = sqrt(1/12)*2*m1
     cens = rnorm(n, m1, m2)   
     status = 1*(tt<=cens)
     y2 = tt*(tt<=cens)+cens*(tt>cens)
     y.o <- log(y2) 
     y.o[which(y.o=="-Inf")] = LOG
     d.o <- 1 * (tt <= cens)      
     xs=x[order(y.o),]
     es=e[order(y.o),]
     d=d.o[order(y.o)]
     yy=sort(y.o)
     w <- numeric(n)
     w[1]=d[1]/n
     for ( i in 2:n )
      {
        tmp = 1
        for ( j in 1: (i-1) )
          tmp = tmp*((n-j)/(n-j+1))^d[j]
     
        w[i]=d[i]/(n-i+1)*tmp
      }

      xbar <- colSums(xs * w)/sum(w)  
      xz <- w * (xs - matrix(rep(xbar, n), ncol=ncol(xs), byrow=TRUE))
      ebar <- colSums(es * w)/sum(w)  
      ez <- w * (es - matrix(rep(ebar, n), ncol=ncol(es), byrow=TRUE))
      Y <- w * yy   

      xx = NULL
           for (i in 1:s)
              {
               temp = xs[,i]*es
               tt.bar <- colSums(temp * w)/sum(w)  
               temp2 <- w * (temp - matrix(rep(tt.bar, n), ncol=ncol(temp), byrow=TRUE))
               xx = cbind(xx,xz[,i],temp2)               
              }
      xx = cbind(ez,xx)

     dat = list(y= Y, x= xz, e= ez, xx=xx, rate=1-sum(status)/n)
     return(dat)    
 
}

