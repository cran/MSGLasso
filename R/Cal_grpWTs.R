Cal_grpWTs <-
function(P, Q, G, R, gmax, PQ.grps)
{
PP=P
QQ=Q
GG=G
RR=R
ggmax=gmax

PQgrps=as.vector(t(PQ.grps))
grpWTs=rep(0, GG*RR)
Gcounts=rep(0,GG)
Rcounts=rep(0,RR)

    junk=.C("Cal_grpWTs",
	    as.integer(PP),
	    as.integer(QQ),
	    as.integer(GG),
	    as.integer(RR),
	    as.integer(ggmax),
	    as.integer(PQgrps),
	    grpWTs=as.double(grpWTs),
	    as.integer(Gcounts),
	    as.integer(Rcounts)
    	    )

grpWTs.result=matrix(junk$grpWTs, nrow=GG, byrow=T)
return(list(grpWTs=grpWTs.result))

}
