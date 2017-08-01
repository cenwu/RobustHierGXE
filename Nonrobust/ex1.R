

          rm(list=ls(all=TRUE))
          library(splines)
          library(glmnet)
          library(MASS)   
          library(mvtnorm)  
          library(quantreg) 
          source("Data1.R")
          source("LShier.r")
          source("LSnhier.r")
          source("LSgroup.r")
          library(mvtnorm)      






          n = 200; s = 200; tol=1e-03; eps=1e-06;          
          reps = 50; nt=10*n;

          rate =rep(0,reps)
          tp.ls = fp.ls = tp.ls1 = fp.ls1 =  tp.ls2 = fp.ls2  =rep(0,reps)
          tp1.ls = fp1.ls = tp1.ls1 = fp1.ls1 = tp1.ls2 = fp1.ls2 = rep(0,reps)
          tp2.ls = fp2.ls = tp2.ls1 = fp2.ls1 = tp2.ls2 = fp2.ls2 = rep(0,reps)
          pred.err.ls = pred.err2.ls = pred.err.ls1 = pred.err2.ls1 = pred.err.ls2 = pred.err2.ls2 =rep(0,reps)


          dat = Data1(n,s)
          y = dat$y; x = dat$x; e = dat$e; xx = dat$xx; 
          # generate an independent testing set
          dat.test = Data1(nt,s)
          y.test = dat.test$y; x.test = dat.test$x; e.test = dat.test$e; xx.test = dat.test$xx;  

          res.ls = LShier(y,x,e,y.test,x.test,e.test,tol,eps,s,xx,xx.test)    
          res.ls1 = LSnhier(y,x,e,y.test,x.test,e.test,tol,eps,s,xx,xx.test)            
          res.ls2 = LSgroup(y,x,e,y.test,x.test,e.test,tol,eps,s,xx,xx.test)   


       
