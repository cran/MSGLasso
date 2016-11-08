#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <Rmath.h>
#include <R.h>

//######################################################//
//######################################################//
//   The Multivariate Sparse Group Lasso      ////////////
//                 and                        ////////////
//   The Mixed Coordinate Descent Algorithm   ////////////
//                                            ////////////
//   Yanming Li                               ////////////
//   2014                                     ////////////
//   liyanmin@umich.edu                       ////////////
//######################################################//
//######################################################//
//
// for the paper "Multivaroate sparse group lasso for the multivaraite multiple linear regression with an arbitrary group structure"
// by Yanming Li, Bin Nan, and Ji Zhu



/////////  I. Accessary functions ///////////////////////////

double Xnorm (int p, int N, int P, double *X_m){
		double Xnorm_p=0;
		for(int n=0;n<N;n++)
		Xnorm_p=Xnorm_p+X_m[n*P+p]*X_m[n*P+p];
		return(Xnorm_p);
}



/////////////
void CalBnorm(int P, int Q, double *Beta, double *Bnorm)
{
	int p,q;
    for(p=0;p<P;p++)
     {
		 Bnorm[p]=0;
         for(q=0;q<Q;q++)
         {
            Bnorm[p]=Bnorm[p]+Beta[p*Q+q]*Beta[p*Q+q];
     }
}
}


/////////////
double CalMBnorm(int P, int Q, double *Beta)
{
    double MBnorm=0;
	int p,q;
    for(p=0;p<P;p++)
     {
         for(q=0;q<Q;q++)
         {
            MBnorm=MBnorm+Beta[p*Q+q]*Beta[p*Q+q];
     }
}
    MBnorm=sqrt(MBnorm);
    return(MBnorm);
}

/////////////
double SoftShrink(double y, double lam)
{
	double temp;
	double result;
	if(y>0)
	 temp=y-lam;
	else
	 temp=-y-lam;
	if(temp<=0)
	 result=0;
	else
	 {
		 result=temp;
		 if(y<0)
		  result=temp*(-1);
	  }
	return(result);
}


//////////////////
void Assign(int P, int Q, double *data_old, double *data_copy)
{
	int p, q;

	for(p=0;p<P;p++)
	 for(q=0;q<Q;q++)
	  data_copy[p*Q+q]=data_old[p*Q+q];
}

/////////////////
void Assignnorm(int P, int Q, double *norm_old, double *norm_copy)
{
	int p, q;

	for(p=0;p<P;p++)
	  norm_copy[p*Q+q]=norm_old[p*Q+q];
}


//////////////////////////////
double Dist(int P, int Q, double *phi_last, double * phi)
{
	int p,q;
	double temp,result;	
      result=0;
	for(p=0;p<P;p++)
	 for(q=0;q<Q;q++)
	  {
		  temp=phi_last[p*Q+q]-phi[p*Q+q];
		  temp=fabs(temp);
		  if(temp>result)
		    result=temp;
      }
    return(result);
}


//////////////////////
void UpdateBeta(int cur_p, int cur_q, int N, int P, int Q, double lambda1, double *lambda_G, double *X_m, double *Xnorm, double *E, double *Beta, double *Beta_old, double *grp_Norm, int G, int R, double *grpWTs, int ggmax, int *PQgrps, int *Pen_L, int *Pen_grp, double *softShrkL2Norm, double *softShrk)
{
	int n, p, q;
	int g, r;
	int i, j;
	double temp, temp1, temp2, tempd, tempsoftShrk;
	p=cur_p;
	q=cur_q;
	int gmax=ggmax;
	double grBnorm_other;
	double grBnorm_other_avg;
	double zeroOtherNormSum;
	int grWeight;
	double grBnorm, grBnorm1;
	double grLambda;


	    temp=0;
	    for(n=0;n<N;n++)
	    temp=temp+E[n*Q+q]*X_m[n*P+p];
	    temp1=temp+Beta_old[p*Q+q]*Xnorm[p];

		if(Pen_L[p*Q+q]==0){zeroOtherNormSum=0;} 
		else{zeroOtherNormSum=lambda1;}
		double nonZeroOtherNormSum=0;

		i=0;		
		while(PQgrps[p*gmax+i]!=999){
			g=PQgrps[p*gmax+i];
		        j=0;
		while(PQgrps[(q+P)*gmax+j]!=999){
			 r=PQgrps[(q+P)*gmax+j];
			 grWeight=grpWTs[g*R+r];
			 if(Pen_grp[g*R+r]!=0){	
			 grBnorm=grp_Norm[g*R+r];
			 if((grBnorm*grBnorm-Beta[p*Q+q]*Beta[p*Q+q]) <0){grBnorm_other=0;}
			 else{grBnorm_other=sqrt(grBnorm*grBnorm-Beta_old[p*Q+q]*Beta_old[p*Q+q]);}

			 grLambda=lambda_G[g*R+r];
			 if ((softShrkL2Norm[g*R+r] <= grWeight*grLambda) || (grBnorm_other< 1e-2)){
			    zeroOtherNormSum=zeroOtherNormSum+lambda_G[g*R+r]*grWeight;}
			 else{   
			    nonZeroOtherNormSum=nonZeroOtherNormSum+lambda_G[g*R+r]*grWeight/grBnorm;}
			}
				 	
		j++;
		}
		i++;
		}
    
	Beta[p*Q+q]=SoftShrink(temp1, N*zeroOtherNormSum)/(Xnorm[p]+ N*nonZeroOtherNormSum);//}
 	tempsoftShrk=SoftShrink(temp1/N, zeroOtherNormSum);

	i=0;
	while(PQgrps[p*gmax+i]!=999){
		g=PQgrps[p*gmax+i];
		j=0;
	while(PQgrps[(q+P)*gmax+j]!=999){		
		r=PQgrps[(q+P)*gmax+j];
		grBnorm=grp_Norm[g*R+r]; 
		if((grBnorm*grBnorm+Beta[p*Q+q]*Beta[p*Q+q]-Beta_old[p*Q+q]*Beta_old[p*Q+q])<0){temp2=0;}
		else{temp2= sqrt(grBnorm*grBnorm+Beta[p*Q+q]*Beta[p*Q+q]-Beta_old[p*Q+q]*Beta_old[p*Q+q]);}
 		grp_Norm[g*R+r]= temp2;

		if((softShrkL2Norm[g*R+r]*softShrkL2Norm[g*R+r]+tempsoftShrk*tempsoftShrk-softShrk[p*Q+q]*softShrk[p*Q+q])<0){tempd=0;}
		else{tempd=sqrt(softShrkL2Norm[g*R+r]*softShrkL2Norm[g*R+r]+tempsoftShrk*tempsoftShrk-softShrk[p*Q+q]*softShrk[p*Q+q]);}
		softShrkL2Norm[g*R+r] =tempd;
	j++;
	}
	i++;
	}


    softShrk[p*Q+q] =tempsoftShrk;

    ////// update residue    
      for(n=0;n<N;n++)
       {
		E[n*Q+q]=E[n*Q+q]+(Beta_old[p*Q+q]-Beta[p*Q+q])*X_m[n*P+p];
	   }


   /////// update BetaLASSO_XY_old    
      Beta_old[p*Q+q]=Beta[p*Q+q];

}



/////////  II Auxilary output functions ///////////////////////////

void Cal_grpWTs(int *PP, int *QQ, int *GG, int *RR, int *ggmax, int *PQgrps, double *grpWTs, int *Gcounts, int *Rcounts)
{
/// P: number of predictor variables
/// Q: number of response variables
/// G: number of covariates groups
/// R: number of response groups
/// gmax: the max number of groups that a signle entry belongs to
/// PQgrps: group indicator matrix of (P+Q) rows and gmax columns
/// grpWTs: group weight matrix of dimension G*R

int P, Q;
int G, R;
int p, q;
int g, r;
int i, j, k;
int gmax;
int pos;

P=*PP;
Q=*QQ;
G=*GG;
R=*RR;
gmax=*ggmax;


for(g=0; g<G; g++){ Gcounts[g]=0;
for(r=0; r<R; r++){ 
Rcounts[r]=0;
 grpWTs[g*R+r]=0;
}
}

for(p=0;p<P;p++)
{
  for(i=0; i<=gmax; i++){
	pos =PQgrps[p*(gmax+1)+i]-0;
	if(pos!=999)
	Gcounts[pos]++;}
}

for(q=0;q<Q;q++)
{
  for(i=0; i<=gmax; i++){
	pos =PQgrps[(P+q)*(gmax+1)+i]-0;
	if(pos!=999)
	Rcounts[pos]++;}
}


for(int i=0; i<G; i++){
  for(int j=0; j<R; j++){
     grpWTs[i*R+j]=sqrt(Gcounts[i]*Rcounts[j]);}
}

}


void Find_PQ_Coord_Grps(int *PP, int *QQ, int *GG, int *RR, int *ggmax, int *GarrStarts, int *GarrEnds, int *RarrStarts, int *RarrEnds, int *PQgrps)
{
/// P: number of predictor variables
/// Q: number of response variables
/// GarrStarts: predictor group start positions
/// GarrEnds: predictor group end positions
/// RarrStarts: response group start positions
/// RarrEnds: responses group End positions
/// gmax: the max number of groups that a signle entry could belong to

int P, Q;
int G, R;
int p, q;
int g, r;
int i, j, k;
int gmax;

P=*PP;
Q=*QQ;
G=*GG;
R=*RR;
gmax=*ggmax+1;


for(i=0; i<(P+Q); i++)
  for(j=0; j<gmax; j++)
     PQgrps[i*gmax+j]=999;

    for(p=0;p<P;p++)
     {
        i=0;
	for(g=0;g<G;g++){
	   if(p>=GarrStarts[g] && p<=GarrEnds[g]){ PQgrps[p*gmax+i] =g; i++;}
	} 
     }

    for(q=0;q<Q;q++)
     {
        j=0;
	for(r=0;r<R;r++){
	   if(q>=RarrStarts[r] && q<=RarrEnds[r]){ PQgrps[(P+q)*gmax+j] =r; j++;}
	} 
     }
}



void Find_GR_Coord_Grps(int *PP, int *QQ, int *GG, int *RR, int *ccmax, int *GarrStarts, int *GarrEnds, int *RarrStarts, int *RarrEnds, int *GRgrps)
{
/// P: number of predictor variables
/// Q: number of response variables
/// GarrStarts: predictor group start positions
/// GarrEnds: predictor group end positions
/// RarrStarts: response group start positions
/// RarrEnds: responses group End positions
/// ccmax: the maximum number of columns or rows a response or predicotr group contains


int P, Q;
int G, R;
int p, q;
int g, r;
int i, j, k;
int cmax;

P=*PP;
Q=*QQ;
G=*GG;
R=*RR;

cmax=*ccmax+1;

for(i=0; i<(G+R); i++)
  for(j=0; j<cmax; j++)
     GRgrps[i*cmax+j]=999;


for(g=0;g<G;g++)
  {
      i=0;
      for(p=0;p<P;p++){
	if(p>=GarrStarts[g] && p<=GarrEnds[g]){ GRgrps[g*cmax+i] =p; i++;}
      } 
  }


for(r=0;r<R;r++)
  {
      j=0;
      for(q=0;q<Q;q++){
	if(q>=RarrStarts[r] && q<=RarrEnds[r]){ GRgrps[(G+r)*cmax+j] =q; j++;}
      } 

   }
}


/////////  Main output functions ///////////////////////////     
/////////  Main function (1) MSGLasso
/////////  without initial input BETA values


void MSGLasso(int *NN, int *PP, int *QQ, int *GG, int *RR,  double *X_m, double *Y_m,  double *grpWTs, int *Pen_L, int *Pen_grp, int *ggmax, int *PQgrps, int *ccmax, int *GRgrps, double *lam1, double *lambda_G, double *grp_Norm, double *Beta_output, int *N_iter, double *RSS, double *E_debug)
{
	/// Variables:
	/// NN: number of observations;
	/// PP: number of predictor variables
	/// QQ: number of response variables
	/// GG: number of predictor groups
	/// RR: number of response groups

	/// X_m: predictor matrix of N by P, input DNA data; each column has mean 0 and sd 1.
	/// Y_m: response matrix of N by Q, input expression data; each column has mean 0.

        /// grpWTs: user specified adaptive group weighting matrix of G by R, for putting different penalization levels on different groups. Can be calculated from the auxilary function Cal_grpWTs.
	/// Pen_L: user specified single-entry level penalization indictor mateix of P by Q. 1 for being penalized and 0 for not.
	/// Pen_grp: user specified group level penalization indictor mateix of G by R. 1 for being penalized and 0 for not.
	/// ggmax: the max number of different groups a single variable belongs to.
	/// PQgrps: the group attribute matrix of (P+Q)*(ggmax+1). See user manual for more details.
	/// ccmax: the max number of variables a single group contains.
	/// GRgrps: the coefficient attribute matrix of (G+R)*(ccmax+1). See user manual for more details.

	/// lam1: lasso panelty parameter;
	/// lambda_G: group penalty parameter;

	/// grp_Norm: a G by R matrix containing the L2 group norms.
	/// Beta_output: matrix of P by Q, estimated coefficients.
	/// N_iter: indicate the total iteration.
	/// RSS: sum of squre error.

	/// E_debug: a intermedia matrix of N by Q, storing the change of the predicted values in each updating step.

	int N, P, Q;
	int G, R;
	int n,p,q;
	int g,r;
	int i,j;
	int n_iter;
	int cur_p, cur_q;
	int gmax;
	int cmax;


        double rss;
	double temp,temp1,temp2;
	double lambda1;
	double *Xnorm;
	double *Beta;
	double *Beta_old, *Beta_last;
	double *E;

	double *softShrk;
	double *softShrkL2Norm;	

	double flag;
	double flag_a;

	double grBnorm, grLambda;
	double grWeight;

	lambda1=*lam1;
	N=*NN;
	P=*PP;
	Q=*QQ;
	G=*GG;
	R=*RR;
	gmax=*ggmax;
	cmax=*ccmax;

   	n_iter=0;

	Xnorm=(double *) malloc (P*sizeof(double));
	Beta=(double *) malloc (P*Q*sizeof(double));
	Beta_old=(double *) malloc (P*Q*sizeof(double));
	Beta_last=(double *) malloc (P*Q*sizeof(double));
	E=(double *) malloc (N*Q*sizeof(double));

	softShrk=(double *) malloc (P*Q*sizeof(double));
	softShrkL2Norm=(double *) malloc (G*R*sizeof(double));


///// initial value
	for(p=0;p<P;p++)
	 {
		 Xnorm[p]=0;
		 for(n=0;n<N;n++)
		 Xnorm[p]=Xnorm[p]+X_m[n*P+p]*X_m[n*P+p];
	 }

//////// Initialize Beta solution
	for(p=0;p<P;p++)
	 for(q=0;q<Q;q++)
	  {
		Beta[p*Q+q]=0;
	  }

//////// Initialize group Norms
	for(g=0;g<G;g++)
	 for(r=0;r<R;r++)
	  {
		grp_Norm[g*R+r]=0; 
	  }

///////// (3) Residue
    for(n=0;n<N;n++)
     for(q=0;q<Q;q++)
      {
		  temp=0;
		  for(p=0;p<P;p++)
		  temp=temp+Beta[p*Q+q]*X_m[n*P+p];
 		  E[n*Q+q]=Y_m[n*Q+q]-temp;
      }


///////// (4) Initialize softShrk

        Assign(P, Q, Beta, Beta_old);

	for(p=0;p<P;p++)
	 for(q=0;q<Q;q++)
	  {
	    temp=0;
	    for(n=0;n<N;n++)
	    temp=temp+E[n*Q+q]*X_m[n*P+p];
	    temp1=temp+Beta_old[p*Q+q]*Xnorm[p];
	    softShrk[p*Q+q]=SoftShrink(temp1/N, Pen_L[p*Q+q]);
	  }


///////// (5) Initialize softShrkL2Norm
		for(g=0;g<G;g++)
		  {
		     for(r=0;r<R;r++)
		       {
				temp=0;
				i=0;
				while(GRgrps[g*cmax+i]!=999){
				  p=GRgrps[g*cmax+i];
		  		  j=0;
		  		  while(GRgrps[(G+r)*cmax+j]!=999){
		      		    q=GRgrps[(G+r)*cmax+j];
				    temp=temp+softShrk[p*Q+q]*softShrk[p*Q+q];
		    		    j++;
		  		  }		  	
		  		  i++;
		                }
			     softShrkL2Norm[g*R+r] =sqrt(temp);
			 }
		   }
		  


        /////////////////////////////////////////////////
        ///////////// (4) begin update
        ////////////////////////////////////////////////

             flag=100;
             while(flag>1e-4 && n_iter<1e+4)
             {

               //////////////// 2) Full loop
               Assign(P, Q, Beta, Beta_last);
               Assign(P, Q, Beta, Beta_old);
               for(p=0;p<P;p++){
		for(q=0;q<Q;q++)               
                {
        		UpdateBeta(p, q, N, P, Q, lambda1, lambda_G, X_m, Xnorm, E, Beta, Beta_old, grp_Norm, G, R, grpWTs, gmax, PQgrps, Pen_L, Pen_grp, softShrkL2Norm, softShrk);

            	} }	

              
              ///////////////  3) Thresholding L2 groups  
            
		for(g=0;g<G;g++)
		  {
		     for(r=0;r<R;r++)
		       {
			  grLambda=lambda_G[g*R+r];
   			  grWeight=grpWTs[g*R+r];
			  if (softShrkL2Norm[g*R+r] <= grWeight*grLambda)
                            {
				i=0;
				while(GRgrps[g*cmax+i]!=999){
				  p=GRgrps[g*cmax+i];
		  		  j=0;
		  		  while(GRgrps[(G+r)*cmax+j]!=999){
		      		    q=GRgrps[(G+r)*cmax+j];
		    		    Beta[p*Q+q]=0;
		    		    j++;
		  		  }		  	
		  		  i++;
		                }
			     grp_Norm[g*R+r] =0;
			    }
		       }
		  }


              n_iter=n_iter+1;
       	      flag=Dist(P, Q, Beta_last, Beta);	
           }//end of all loop while(flag>1e-4 && n_iter<1e+4)

//////////////////////////


//////////// calculate Residue
    rss=0;  
    for(n=0;n<N;n++) 
     for(q=0;q<Q;q++)
	   {
		  rss=rss+E[n*Q+q]*E[n*Q+q];
	   }
    Assign(N, Q, E, E_debug);

//////////// return Beta
Assign(P, Q, Beta, Beta_output);
*N_iter=n_iter;
*RSS=rss;

//////// free allocated variables
free(Xnorm);
free(Beta);
free(Beta_old);
free(Beta_last);
free(E);
}///end MSGLasso function 


/////////  Main output functions ///////////////////////////     
/////////  Main function (2) MSGLasso_Ini
/////////  with initial input BETA values

void MSGLasso_Ini(int *NN, int *PP, int *QQ, int *GG, int *RR, double *X_m, double *Y_m,  double *grpWTs, int *Pen_L, int *Pen_grp, int *ggmax, int *PQgrps, int *ccmax, int *GRgrps, double *lam1, double *lambda_G, double *Beta_initial, double *grp_Norm, double *Beta_output, int *N_iter, double *RSS, double *E_debug)
{
	/// Variables:
	/// NN: number of observations;
	/// PP: number of predictor variables
	/// QQ: number of response variables
	/// GG: number of predictor groups
	/// RR: number of response groups

	/// X_m: predictor matrix of N by P, input DNA data; each column has mean 0 and sd 1.
	/// Y_m: response matrix of N by Q, input expression data; each column has mean 0.

        /// grpWTs: user specified adaptive group weighting matrix of G by R, for putting different penalization levels on different groups. Can be calculated from the auxilary function Cal_grpWTs.
	/// Pen_L: user specified single-entry level penalization indictor mateix of P by Q. 1 for being penalized and 0 for not.
	/// Pen_grp: user specified group level penalization indictor mateix of G by R. 1 for being penalized and 0 for not.
	/// ggmax: the max number of different groups a single variable belongs to.
	/// PQgrps: the group attribute matrix of (P+Q)*(ggmax+1). See user manual for more details.
	/// ccmax: the max number of variables a single group contains.
	/// GRgrps: the coefficient attribute matrix of (G+R)*(ccmax+1). See user manual for more details.

	/// lam1: lasso panelty parameter;
	/// lambda_G: group penalty parameter;

	/// grp_Norm: a G by R matrix containing the L2 group norms.
	/// Beta_output: matrix of P by Q, estimated coefficients.
	/// N_iter: indicate the total iteration.
	/// RSS: sum of squre error.

	/// E_debug: a intermedia matrix of N by Q, storing the change of the predicted values in each updating step.

	int N, P, Q;
	int G, R;
	int n,p,q;
	int g,r;
	int i,j;
	int n_iter;
	int cur_p, cur_q;
	int gmax;
	int cmax;
 

    	double rss;
	double temp,temp1,temp2;
	double lambda1;
	double *Xnorm;
	double *Beta;
	double *Beta_old, *Beta_last;
	double *E;

	double *softShrk;
	double *softShrkL2Norm;	

	double flag;
	double flag_a;

	double grBnorm, grLambda;
	double grWeight;
	double grSize;

	lambda1=*lam1;
	N=*NN;
	P=*PP;
	Q=*QQ;
	G=*GG;
	R=*RR;
	gmax=*ggmax;
 	cmax=*ccmax;

   	n_iter=0;

	Xnorm=(double *) malloc (P*sizeof(double));

	Beta=(double *) malloc (P*Q*sizeof(double));
	Beta_old=(double *) malloc (P*Q*sizeof(double));
	Beta_last=(double *) malloc (P*Q*sizeof(double));
	E=(double *) malloc (N*Q*sizeof(double));

	softShrk=(double *) malloc (P*Q*sizeof(double));
	softShrkL2Norm=(double *) malloc (G*R*sizeof(double));

///// initial value
	for(p=0;p<P;p++)
	 {
		 Xnorm[p]=0;
		 for(n=0;n<N;n++)
		 Xnorm[p]=Xnorm[p]+X_m[n*P+p]*X_m[n*P+p];
	 }


//////// Initialize Beta solution
       for(p=0;p<P;p++)
       for(q=0;q<Q;q++)
         {
	   Beta[p*Q+q]=Beta_initial[p*Q+q];
         }

///////// (3) Residue
    for(n=0;n<N;n++)
     for(q=0;q<Q;q++)
      {
		  temp=0;
		  for(p=0;p<P;p++)
		  temp=temp+Beta[p*Q+q]*X_m[n*P+p];
 		  E[n*Q+q]=Y_m[n*Q+q]-temp;
      }


///////// (4) Initialize softShrk

        Assign(P, Q, Beta, Beta_old);

	for(p=0;p<P;p++)
	 for(q=0;q<Q;q++)
	  {
	    temp=0;
	    for(n=0;n<N;n++)
	    temp=temp+E[n*Q+q]*X_m[n*P+p];
	    temp1=temp+Beta_old[p*Q+q]*Xnorm[p];
	    softShrk[p*Q+q]=SoftShrink(temp1/N, Pen_L[p*Q+q]);
	  }


///////// (5) Initialize softShrkL2Norm
		for(g=0;g<G;g++)
		  {
		     for(r=0;r<R;r++)
		       {
				temp=0;
				i=0;
				while(GRgrps[g*cmax+i]!=999){
				  p=GRgrps[g*cmax+i];
		  		  j=0;
		  		  while(GRgrps[(G+r)*cmax+j]!=999){
		      		    q=GRgrps[(G+r)*cmax+j];
				    temp=temp+softShrk[p*Q+q]*softShrk[p*Q+q];
		    		    j++;
		  		  }		  	
		  		  i++;
		                }
			     softShrkL2Norm[g*R+r] =sqrt(temp);
			 }
		    }
		  


        /////////////////////////////////////////////////
        ///////////// (4) begin update
        ////////////////////////////////////////////////

             flag=100;
             while(flag>1e-4 && n_iter<1e+4)
             {

               //////////////// 2) Full loop
               Assign(P, Q, Beta, Beta_last);
               Assign(P, Q, Beta, Beta_old);
               for(p=0;p<P;p++){
		for(q=0;q<Q;q++)               
                {
        		UpdateBeta(p, q, N, P, Q, lambda1, lambda_G, X_m, Xnorm, E, Beta, Beta_old, grp_Norm, G, R, grpWTs, gmax, PQgrps, Pen_L, Pen_grp, softShrkL2Norm, softShrk);
	
            	} }


              ///////////////  3) Thresholding L2 groups 
             
		for(g=0;g<G;g++)
		  {
		     for(r=0;r<R;r++)
		       {
			  grLambda=lambda_G[g*R+r];
   			  grWeight=grpWTs[g*R+r];
			  if (softShrkL2Norm[g*R+r] <= grWeight*grLambda)
                           {
				i=0;
				while(GRgrps[g*cmax+i]!=999){
				  p=GRgrps[g*cmax+i];
		  		  j=0;
		  		  while(GRgrps[(G+r)*cmax+j]!=999){
		      		    q=GRgrps[(G+r)*cmax+j];
		    		    Beta[p*Q+q] =0;
		    		    j++;
		  		  }		  	
		  		  i++;
		                }
				grp_Norm[g*R+r] =0;
			    }
		       }
		  }



	n_iter=n_iter+1;	
       	flag=Dist(P, Q, Beta_last, Beta);	
            }//end of all loop while(flag>1e-4 && n_iter<1e+4)


//////////////////////////


//////////// calculate Residue
    rss=0;  
    for(n=0;n<N;n++) 
     for(q=0;q<Q;q++)
	   {
		  rss=rss+E[n*Q+q]*E[n*Q+q];
	   }
    Assign(N, Q, E, E_debug);

//////////// return Beta
Assign(P, Q, Beta, Beta_output);
*N_iter=n_iter;
*RSS=rss;

//////// free allocated variables
free(Xnorm);
free(Beta);
free(Beta_old);
free(Beta_last);
free(E);
}///end MSGLasso function 



