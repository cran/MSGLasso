FindingPQGrps <-
function(P, Q, G, R, gmax, G.Starts, G.Ends, R.Starts, R.Ends)
{
PP=P
QQ=Q
GG=G
RR=R

ggmax=gmax
GarrStarts=G.Starts
GarrEnds=G.Ends
RarrStarts=R.Starts
RarrEnds=R.Ends

PQgrps <- matrix(0, nrow=(PP+QQ), ncol=ggmax+1)

   junk=.C("Find_PQ_Coord_Grps",
          as.integer(PP),
          as.integer(QQ),
          as.integer(GG),
	  as.integer(RR),

	  as.integer(ggmax),
	  as.integer(GarrStarts),
	  as.integer(GarrEnds),

	  as.integer(RarrStarts),
	  as.integer(RarrEnds),
	  PQgrps.out=as.integer(PQgrps)
          )  

PQgrps.result=matrix(junk$PQgrps.out, nrow=(PP+QQ), byrow=T)
return(list(PQgrps=PQgrps.result))

}
