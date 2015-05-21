MSGLasso <-
function(X.m, Y.m, grp.WTs, Pen.L, Pen.G, PQ.grps, GR.grps, grp_Norm0, lam1, lam.G, Beta0=NULL)
{
##########################################################

  N=nrow(X.m)
  P=ncol(X.m)
  if(!is.null(ncol(Y.m))){
   Q=ncol(Y.m)} else {
   Q=1
  }
	
  G=nrow(grp.WTs)
  if(!is.null(ncol(grp.WTs))){
   R=ncol(grp.WTs)} else {
   R=1 
  }

  X.v=as.vector(t(X.m))
  if(!is.null(ncol(Y.m))){
   Y.v=as.vector(t(Y.m))} else {
   Y.v=as.vector(Y.m)
  }

  if(!is.null(ncol(grp.WTs))){
   grpWTs=as.vector(t(grp.WTs))} else{
   grpWTs=as.vector(grp.WTs)
  }

  if(!is.null(ncol(Pen_L))){
   Pen_L =as.vector(t(Pen.L))} else{
   Pen_L =as.vector(Pen.L)
  }
  if(!is.null(ncol(Pen_G))){
   Pen_G =as.vector(t(Pen.G))} else{
   Pen_G =as.vector(Pen.G)   
  }

  gmax=ncol(PQ.grps)
  if(is.null(gmax) | gmax==1){
    err <- sprintf("PQ.grps should be a matrix of (P+Q) rows and (gmax+1) columns.")
    stop(err)
  } else {
    PQgrps=as.vector(t(PQ.grps))
  }

  cmax=ncol(GR.grps)
  if(is.null(cmax) | cmax==1){
    err <- sprintf("GR.grps should be a matrix of (G+R) rows and (cmax+1) columns.")
    stop(err)
  } else {
    GRgrps=as.vector(t(GR.grps))
  }

  
  lambda1=lam1
  lambdaG=as.vector(t(lam.G))	

  Beta.m=matrix(0, nrow=P, ncol=Q)
  Beta.v=as.vector(Beta.m)
  grp_Norm.v=as.vector(t(grp_Norm0))  

  iter.count=0
  RSS=0
  Edebug=rep(0, N*Q)


  ###### begin estimation
  if(is.null(Beta0))
  {
   junk=.C("MSGLasso",
          as.integer(N),
          as.integer(P),
          as.integer(Q),
          as.integer(G),
          as.integer(R),

          as.double(X.v),
          as.double(Y.v),

	  as.double(grpWTs),

	  as.integer(Pen_L),
	  as.integer(Pen_G),

	  as.integer(gmax),
	  as.integer(PQgrps),
	  as.integer(cmax),
	  as.integer(GRgrps),

          as.double(lambda1),
	  as.double(lambdaG),

	  grp_Norm=as.double(grp_Norm.v),
          Beta.out=as.double(Beta.v),
          n.iter=as.integer(iter.count),
          RSS=as.double(RSS),
          Edebug=as.double(Edebug)
          )  
  } else {
       Beta.ini.v=as.vector(t(Beta0))
   junk=.C("MSGLasso_Ini",
          as.integer(N),
          as.integer(P),
          as.integer(Q),
          as.integer(G),
          as.integer(R),

          as.double(X.v),
          as.double(Y.v),

	  as.double(grpWTs),

	  as.integer(Pen_L),
	  as.integer(Pen_G),

	  as.integer(gmax),
	  as.integer(PQgrps),
	  as.integer(cmax),
	  as.integer(GRgrps),

          as.double(lambda1),
	  as.double(lambdaG),

 	  Beta.ini=as.double(Beta.ini.v),
	  grp_Norm=as.double(grp_Norm.v),
          Beta.out=as.double(Beta.v),
          n.iter=as.integer(iter.count),
          RSS=as.double(RSS),
          Edebug=as.double(Edebug)

          )  
  }


  Beta.result=matrix(junk$Beta.out, nrow=P, byrow=T)
  grp_Norm.result=matrix(junk$grp_Norm, nrow=G, byrow=T)
  E=matrix(junk$Edebug, nrow=N, ncol=Q, byrow=T)
  rss.v=apply(E^2,1,sum)
  rss=sum(rss.v)
  iter.count=junk$n.iter

  return(list(Beta=Beta.result, grpNorm=grp_Norm.result, rss.v=rss.v, rss=rss, iter=iter.count))
}
