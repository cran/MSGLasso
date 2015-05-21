FindingGRGrps <-
function(P, Q, G, R, cmax, G.Starts, G.Ends, R.Starts, R.Ends)
{
PP=P
QQ=Q
GG=G
RR=R

ccmax=cmax
GarrStarts=G.Starts
GarrEnds=G.Ends
RarrStarts=R.Starts
RarrEnds=R.Ends

GRgrps <- matrix(0, nrow=(GG+RR), ncol=(ccmax+1))

   junk=.C("Find_GR_Coord_Grps",
          as.integer(PP),
          as.integer(QQ),
          as.integer(GG),
	  as.integer(RR),

	  as.integer(ccmax),
	  as.integer(GarrStarts),
	  as.integer(GarrEnds),

	  as.integer(RarrStarts),
	  as.integer(RarrEnds),
	  GRgrps.out=as.integer(GRgrps)
          )  

GRgrps.result=matrix(junk$GRgrps.out, nrow=(GG+RR), byrow=T)
return(list(GRgrps=GRgrps.result))

}
