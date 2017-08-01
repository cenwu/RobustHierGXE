

          rm(list=ls(all=TRUE))
          library(splines)
          library(glmnet)
          library(MASS)   
          library(mvtnorm)  
          library(quantreg) 
          source("Data1.R")
          source("LADhier.r")
          source("LADnhier.r")
          source("LADgroup.r")
          library(mvtnorm)      
          library(splines)
          library(glmnet)
          library(quantreg) 


          n = 200; s = 200; tol=1e-03; eps=1e-06;          
          reps = 50; nt=10*n;
          rate =rep(0,reps)
          tp.lad = fp.lad = tp.lad1 = fp.lad1 =  tp.lad2 = fp.lad2  =rep(0,reps)
          tp1.lad = fp1.lad = tp1.lad1 = fp1.lad1 = tp1.lad2 = fp1.lad2 = rep(0,reps)
          tp2.lad = fp2.lad = tp2.lad1 = fp2.lad1 = tp2.lad2 = fp2.lad2 = rep(0,reps)
          pred.err.lad = pred.err2.lad = pred.err.lad1 = pred.err2.lad1 = pred.err.lad2 = pred.err2.lad2 =rep(0,reps)


          dat = Data1(n,s)
          y = dat$y; x = dat$x; e = dat$e; xx = dat$xx; 
          # generate an independent testing set
          dat.test = Data1(nt,s)
          y.test = dat.test$y; x.test = dat.test$x; e.test = dat.test$e; xx.test = dat.test$xx;  

          res.lad = LADhier(y,x,e,y.test,x.test,e.test,tol,eps,s,xx,xx.test)              
          res.lad1 = LADnhier(y,x,e,y.test,x.test,e.test,tol,eps,s,xx,xx.test)          
          res.lad2 = LADgroup(y,x,e,y.test,x.test,e.test,tol,eps,s,xx,xx.test)   




