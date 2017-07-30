

solve.new <- function(mat)
{
  k = nrow(mat)
  tiny = matrix(0,k,k)
  diag(tiny) = 0.1*rep(1,k)
  p = svd(mat+tiny)
  f = p$d
  f[f <= 0.001] = 0.001
  m = matrix(0,k,k)
  diag(m) = f^(-1)
  inv = p$v%*%m%*%t(p$u)
return(inv)
}

galasso<-function(norm,lambda)
{

       if (norm <= lambda)
         res = 0
       else 
         res = norm - lambda
       return(res)
}   


LSgroup <- function(y,x,e,y.test,x.test,e.test,tol,eps,s,xx,xx.test)  
{     

           lambda1 =   seq(180,600, by = 150)/(1*10^14) ; nlambda1 = length(lambda1);
           lambda2 =   lambda1 ; nlambda2 = 1
           n = length(y); nt=10*n; tau = 0.5; q = 10
           design = NULL
           res = matrix(0,n,1)    
           pred.error2 = sum((y.test-mean(y.test))^2)
           beta0 = matrix(0,dim(xx)[2],1) 
           lasso.cv <- cv.glmnet(xx,y,alpha=1,nfolds=5)   
           lasso.fit1 <- glmnet(xx,y,family="gaussian",alpha=1,nlambda=100,intercept=FALSE)  
           beta0 <- as.matrix(as.vector(predict(lasso.fit1, s=lasso.cv$lambda.min, type="coefficients"))[-1]) 
           res0 = y - as.matrix(xx) %*% as.matrix(beta0)
           y.new = y - as.matrix(xx[,((q+1):dim(xx)[2])])%*%as.matrix(beta0[(q+1):dim(xx)[2]])
           lm = lm(y.new~xx[,1:q]-1)
           coef.main = lm$coefficients
           beta0[1:q] = coef.main 
           p1 = dim(beta0)[1]
           p2 = dim(beta0)[2]
           coef = array(0, dim=c(p1,p2,nlambda1,nlambda2))
           cv.test = array(0, dim=c(nlambda1,nlambda2))
           cv.train = array(0, dim=c(nlambda1,nlambda2))
           w.g1 = matrix(0,s,1)
           beta = beta0;
           for(m in 1:s) 
               {
                  sub = (q+(q+1)*(m-1)+1):(q+(q+1)*m)
                  x.sub = xx[,sub]
                  bm = beta[sub] 
                  norm.sub = sqrt(mean((x.sub%*%bm)^2))
                  w.g1[m] = (max(norm.sub,eps/100))^(-2)  
               }
          w.g2 =  (abs(beta)+eps/100)^(-1)  

          for (i in 1:nlambda1)
             {       
 
              for (j in 1:nlambda2)
                {      
                   lam1 = lambda1[i]  
                   lam2 = 0  
                   if (i==1) {  beta = beta0 }  
                   beta = beta0                  
                   res = res0                  
                   count = array(0, dim=c(n,nlambda1))
                   diff.out.old = 2 
                   diff.in = diff.out = 1 
                   step.in = step.out = 0
                   while(diff.out>=10^{-4} & step.out<=15 ) 
                      {
                            beta.out = beta
                                 for(m in 1:s) 
                                     {
                                            sub = (q+(q+1)*(m-1)+1):(q+(q+1)*m)
                                            y.sub = y-xx[,-sub]%*%beta[-sub]
                                            x.sub = xx[,sub]
                                            bm = beta[sub]                                      
                                            w.g3 = (abs(bm)+eps/100)^(-1) 
                                            um = (t(x.sub)%*%y.sub)/n 
                                            vm = as.matrix(sign(um)*(abs(um)-lam2*w.g3)*I(abs(um)>lam2*w.g3))
                                            sm = solve.new(t(x.sub)%*%x.sub/n)%*%vm
                                            norm.sub = sqrt(mean((x.sub%*%bm)^2))
                                            norm.sm = sqrt(mean((x.sub%*%sm)^2))
                                            if (sum(abs(bm))>0)
                                              {
                                                if (sum(abs(vm))>0) {  beta[sub] = galasso(norm.sm,lam1*w.g1[m])/norm.sm*sm }
                                                if (sum(abs(vm))==0) { beta[sub] = vm }
                                              }                                    

                                     }      
                          beta[abs(beta)<0.00001] <- 0                     
                          diff.out.old = diff.out;                         
                          diff.out = mean((beta-beta.out)^2); 
                          step.out = step.out+1                        
                          cat("step.out = ", step.out, "lam1 = ", lam1, "\n")
                          cat("diff.out = ", diff.out,  "\n")
                          #print(j)
                          y.new = y - as.matrix(xx[,((q+1):dim(xx)[2])])%*%as.matrix(beta[(q+1):dim(xx)[2]])
                          lm = lm(y.new~xx[,1:q]-1)
                          coef.main = lm$coefficients
                          beta[1:q] = coef.main                       
                       }   
                   res.test = y.test - as.matrix(xx.test) %*% as.matrix(beta)                      
                   cv.test[i,j] =  sum(res.test^2) 
                   lm.test = lm(y.test~xx.test[,(beta!=0)]-1)
                   cv.test[i,j] = sum(lm.test$residuals^2)
                   res = y - as.matrix(xx) %*% as.matrix(beta)     
                   cv.train[i,j] =  sum(res.test^2)    
                   coef[,,i,j] = beta   
             } 
           } 
          cv.train = cv.train/n
          cv.test = cv.test/nt
          cv = cv.test
          i1.cv = which(cv == min(cv), arr.ind = TRUE)[1,1]
          i2.cv = 1; #which(cv == min(cv), arr.ind = TRUE)[1,2]
          flag.lambda1 = lambda1[i1.cv]          
          flag.lambda2 = lambda2[i2.cv]     
          pred.error = cv[i1.cv,i2.cv]*nt
          beta.p = coef[,,i1.cv,i2.cv]  
          ind.m = seq(11,110,by=11) 
          ind.i = c(12,14,16,23,25,27,34,36,38,45,47,49,56,58,60) 
          ind.m1 = seq(121,2200,by=11)
          ind.i1 = seq(1,2210,by=1)[-c(1:q,ind.m,ind.i,ind.m1)]
          TP1 = sum((which(beta.p!=0) %in% ind.m)) 
          TP2 = sum((which(beta.p!=0) %in% ind.i)) 
          TP = TP1 + TP2
          FP1 = sum((which(beta.p!=0) %in% ind.m1))
          FP2 = sum((which(beta.p!=0) %in% ind.i1))
          FP = FP1 + FP2

                     
tun = list(beta.p=beta.p,pred.error=pred.error,pred.error2=pred.error2, TP=TP, FP = FP, TP1=TP1, FP1=FP1,TP2=TP2, FP2=FP2, cv=cv.test,cv2=cv.train)
return(tun)
}  




