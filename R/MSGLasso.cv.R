MSGLasso.cv <-
function(X, Y, grpWTs, Pen.L, Pen.G, PQgrps, GRgrps, lam1.v, lamG.v, fold=10, seed=1, Beta.ini=NULL, grp_Norm=NULL)
{

  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  G=nrow(Pen.G)
  R=ncol(Pen.G)

 ###### sort penalty parameter from large to small
 k1=length(lam1.v) 
 lambda1.v=sort(lam1.v)[k1:1]

 k3=length(lamG.v)
 lambda3.v=sort(lamG.v)[k3:1]
    
 ############################### set seed and generate cross validation seperations
 set.seed(seed)

 index.cv<-NULL
 ran.order=sample(1:n, n, replace=F)
 f.n=floor(n/fold)
 for (f in 1:(fold-1))
    {  
     index.cv[[f]]<-ran.order[(f-1)*f.n+1:f.n] 
    } 
 index.cv[[fold]]<-ran.order[((fold-1)*f.n+1):n]
 
 ################################ begin cross validation

 rss.cv<-array(0,dim=c(k1,k3))

 lams.c <- NULL
 if(is.null(grp_Norm)){
  grp_Norm0 <- matrix(rep(0, G*R), nrow=G, byrow=T)}
 else{
  grp_Norm0 <- grp_Norm
 }
 

 for(f in 1:fold)
  {  print(paste("fold=",f))   
     index.cur.cv<-index.cv[[f]]   
     X.m<-X[-(index.cur.cv),]   
     Y.m<-Y[-(index.cur.cv),]   
     X.t<-X[index.cur.cv,]   
     Y.t<-Y[index.cur.cv,]   

     rss.cv.c<-array(0,dim=c(k1,k3))
   
     Beta.old <-Beta.ini

     for(j in 1:k3)
     {      
       
       if(j>1) Beta.old <- Beta.lasti
       
       for(i in 1:k1)       
       { 

          print(paste(i,j))               
          cur.lam1=lambda1.v[i]          
          cur.lam3=lambda3.v[j]
	  lam1=cur.lam1
	  lam_G=matrix(rep(cur.lam3, G*R),G,R,byrow=T)

          temp=MSGLasso(X.m, Y.m, grpWTs, Pen.L, Pen.G, PQgrps, GRgrps, grp_Norm0, lam1, lam_G, Beta0=Beta.old)   
           
          lams.c[[(j-1)*k1+i]]<-list(lam1=cur.lam1, lam3=cur.lam3) #, phi=temp$Beta#
          Beta.old <-temp$Beta
	  grp_Norm0 <- temp$grpNorm
          
          rss.cv.c[i,j]<-RSS.CV(X.t,Y.t,temp$Beta)
         
          if(i==1) Beta.lasti=temp$Beta 
        } ### end of i iter     
      }### end of j iter

       rss.cv<-rss.cv+rss.cv.c
    
   }###end fold loop
   
   rss.cv<-rss.cv/fold

   ################################## return CV results
   result=list(rss.cv=rss.cv, lams.c=lams.c)
   return(result)
}


####################################################

RSS.CV<-function(X.m,Y.m, phi.est)
{    
    Y.est<-X.m%*%phi.est       
    resi<-Y.m-Y.est
    rss<-sum(as.vector(resi)^2)
    return(rss)
}


