#include "udf.h"
#include <stdio.h>
#include <math.h>
#include <cxiface.h>
#include "surf.h"
#include "para.h"
#include "dpm.h"
real points[6][22311];
typedef struct
{
  Thread *t; 
  cell_t c;
  real xc[ND_ND]; 
  int niptype; 
  real interplow[ND_ND];
  real interphigh[ND_ND];
  real lenghtofinterface;
  real vof;
  real source;
} domainmatrix;
domainmatrix mymatrix[67][333];
real D=0.00000000027, Di=0.000000001;
real b=0.000663;  real gaminfy=0.00000291;

DEFINE_EXECUTE_AT_END(DomainCreater)/*Domain creater produces a matrix just like the geometry and puts in every cell the cell cell_t, thread and the cell centroid*/
{
  Domain *d;  Thread *t;  cell_t c;  real xc[ND_ND]; 
  d=Get_Domain(2);
  int i=332;int j=0; 
  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
	{
	  mymatrix[j][i].t=t;	  mymatrix[j][i].c=c;	  C_CENTROID(xc,c,t);	  mymatrix[j][i].xc[0]=xc[0];	  mymatrix[j][i].xc[1]=xc[1]; mymatrix[j][i].vof=C_VOF(c,t);
	 
	  i=i-1;
	  if(i==-1)
	    {
	      i=332;
	      j=j+1;
	    }
	}
      end_c_loop(c,t)
	}
}
DEFINE_EXECUTE_AT_END(Domainchecker)/*Domain creater produces a matrix just like the geometry and puts in every cell the cell cell_t, thread and the cell centroid*/
{
  Domain *d;  Thread *t;  cell_t c;  real xc[ND_ND], xc1[ND_ND]; int count=0;
  d=Get_Domain(2);
  int i=332,j=0; 
  FILE *fp; fp=fopen("Domain check.txt","a");
  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
	{
	  C_CENTROID(xc,c,t);
	  /* fprintf(fp,"current cell : ( %g,  %g)             ", xc[0], xc[1]);*/
	   if (i>=1 && i<=331 && j<=65 && j>=1)
	    {	 
	      	fprintf(fp,"vof : %g\n",mymatrix[j][i].vof);
		if (mymatrix[j][i].vof==0)
		{
		  if (mymatrix[j+1][i].vof==1 || mymatrix[j-1][i].vof==1 || mymatrix[j][i+1].vof==1 || mymatrix[j][i-1].vof==1)
		    { 
		      /* fprintf(fp,"Interface passesm thru a face.");*/ count++;
		    }
		}
		}
	   /* fprintf(fp,"\n");*/
	  
	  /*	  if (i==0 && j==0){  fprintf(fp,"corner 4 detected\n");}
	  if (i==0 && j==66){ fprintf(fp,"corner 1 detected\n");}
	  if (i==332 && j==0){ fprintf(fp,"corner 3 detected\n");}
	  if (i==332 && j==66){ fprintf(fp,"corner 2 detected\n");}*/
	  i=i-1;
	  if(i==-1)
	    {
	      i=332;
	      j=j+1;
	    }
	}
      end_c_loop(c,t)
	}
  if (count!=0)/*change*/
    {
      fprintf(fp,"Number of interface cells misbehaving is : %d\n\n",count);
    }
  fclose(fp);
}

DEFINE_EXECUTE_AT_END(InterfaceManager)/*can this be patented*/
{ 
  int i=332,j=0; 
  Domain *d;  Thread *t;  cell_t c;  real xc[ND_ND];
  real x1[ND_ND],    x2[ND_ND],     x3[ND_ND],     x4[ND_ND];
  real vof1,vof2,vof3,vof4;  d=Get_Domain(2);         
  real nc=0;/*nc is the number of neibhours of the current cell. should be made 0 in begginning of ex=very new cell.*/
  int nip=0;/*nip is the number of interface points in the current cell. It should be made zero for every new cell.*/
  
  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
	{
	  vof1=0; vof2=0; vof3=0; vof4=0;	  nc=0;	  nip=0;	  mymatrix[j][i].niptype=0;
	  mymatrix[j][i].lenghtofinterface=0;/*change*/
	  C_CENTROID(xc,c,t);
	  x1[0] =xc[0]-0.0003003/2;	  x1[1] =xc[1]+0.0003003/2;	
	  x2[0] =xc[0]+0.0003003/2;	  x2[1] =xc[1]+0.0003003/2;	 
	  x3[0] =xc[0]+0.0003003/2;	  x3[1] =xc[1]-0.0003003/2;	 
	  x4[0] =xc[0]-0.0003003/2;       x4[1] =xc[1]-0.0003003/2;
	  if (i>=1 && i<=331 && j<=65 && j>=1)
	    {  
	      vof1=mymatrix[j+1][i-1].vof + mymatrix[j][i].vof + mymatrix[j+1][i].vof + mymatrix[j][i-1].vof;
	      vof2=mymatrix[j+1][i+1].vof + mymatrix[j][i].vof + mymatrix[j+1][i].vof + mymatrix[j][i+1].vof;
	      vof3=mymatrix[j-1][i+1].vof + mymatrix[j][i].vof + mymatrix[j-1][i].vof + mymatrix[j][i+1].vof;
	      vof4=mymatrix[j-1][i-1].vof + mymatrix[j][i].vof + mymatrix[j-1][i].vof + mymatrix[j][i-1].vof;
	      vof1=vof1/4; vof2=vof2/4; vof3=vof3/4; vof4=vof4/4;
	    }
	  else if(i==0 && j>=1 && j<=65)
	    {  
	      vof1= mymatrix[j][i].vof + mymatrix[j+1][i].vof;
	      vof2=mymatrix[j+1][i+1].vof + mymatrix[j][i].vof + mymatrix[j+1][i].vof + mymatrix[j][i+1].vof;
	      vof3=mymatrix[j-1][i+1].vof + mymatrix[j][i].vof + mymatrix[j-1][i].vof + mymatrix[j][i+1].vof;
	      vof4= mymatrix[j][i].vof + mymatrix[j-1][i].vof;
	      vof1=vof1/2; vof2=vof2/4; vof3=vof3/4; vof4=vof4/2;
	    }
	  else if(i==332 && j>=1 && j<=65)
	    {	 
	      vof1=mymatrix[j+1][i-1].vof + mymatrix[j][i].vof + mymatrix[j+1][i].vof + mymatrix[j][i-1].vof;
	      vof2= mymatrix[j][i].vof + mymatrix[j+1][i].vof;
	      vof3= mymatrix[j][i].vof + mymatrix[j-1][i].vof;
	      vof4=mymatrix[j-1][i-1].vof + mymatrix[j][i].vof + mymatrix[j-1][i].vof + mymatrix[j][i-1].vof;
	      vof1=vof1/4;vof2=vof2/2;vof3=vof3/2; vof4=vof4/4;
	    }
	  else if (j==66 && i<=331 && i>=1)
	    {	
	      vof1= mymatrix[j][i].vof + mymatrix[j][i-1].vof;
	      vof2= mymatrix[j][i].vof + mymatrix[j][i+1].vof;
	      vof3=mymatrix[j-1][i+1].vof + mymatrix[j][i].vof + mymatrix[j-1][i].vof + mymatrix[j][i+1].vof;
	      vof4=mymatrix[j-1][i-1].vof + mymatrix[j][i].vof + mymatrix[j-1][i].vof + mymatrix[j][i-1].vof;
	      vof1=vof1/2; vof2=vof2/2; vof3=vof3/4; vof4=vof4/4;
	      
	    }
	  else if (j==0 && i<=331 && i>=1)
	    {	 
	      vof1=mymatrix[j+1][i-1].vof + mymatrix[j][i].vof + mymatrix[j+1][i].vof + mymatrix[j][i-1].vof;
	      vof2=mymatrix[j+1][i+1].vof + mymatrix[j][i].vof + mymatrix[j+1][i].vof + mymatrix[j][i+1].vof;
	      vof3= mymatrix[j][i].vof + mymatrix[j][i+1].vof;
	      vof4= mymatrix[j][i].vof + mymatrix[j][i-1].vof;
	      vof1=vof1/4; vof2=vof2/4; vof3=vof3/2; vof4=vof4/2;
	     
	    }
	  else if(j==0 && i==0)
	    {	 
	      vof1= mymatrix[j][i].vof + mymatrix[j+1][i].vof;
	      vof2= mymatrix[j+1][i+1].vof + mymatrix[j][i].vof + mymatrix[j+1][i].vof + mymatrix[j][i+1].vof;
	      vof3= mymatrix[j][i].vof + mymatrix[j][i+1].vof;
	      vof4= mymatrix[j][i].vof;
	      vof1=vof1/2; vof2=vof2/4; vof3=vof3/2; vof4=vof4/1;

	    }
	  else if(j==0 && i==332)
	    {	  
	      vof1=mymatrix[j+1][i-1].vof + mymatrix[j][i].vof + mymatrix[j+1][i].vof + mymatrix[j][i-1].vof;
	      vof2=mymatrix[j][i].vof + mymatrix[j+1][i].vof;
	      vof3= mymatrix[j][i].vof;
	      vof4= mymatrix[j][i].vof + mymatrix[j][i-1].vof;
	      vof1=vof1/4; vof2=vof2/2; vof3=vof3/1; vof4=vof4/2;
	     
	    }
	  else if(j==66 && i==0)
	    {	 
	      vof1= mymatrix[j][i].vof;
	      vof2= mymatrix[j][i].vof + mymatrix[j][i+1].vof;
	      vof3=mymatrix[j-1][i+1].vof + mymatrix[j][i].vof + mymatrix[j-1][i].vof + mymatrix[j][i+1].vof;
	      vof4= mymatrix[j][i].vof + mymatrix[j-1][i].vof;
	      vof1=vof1/1; vof2=vof2/2; vof3=vof3/4; vof4=vof4/2;
	    }
	  else if (j==66 && i==332)
	    {  
	      vof1= mymatrix[j][i].vof + mymatrix[j][i-1].vof;
	      vof2= mymatrix[j][i].vof;
	      vof3= mymatrix[j][i].vof + mymatrix[j-1][i].vof;
	      vof4=mymatrix[j-1][i-1].vof + mymatrix[j][i].vof + mymatrix[j-1][i].vof + mymatrix[j][i-1].vof;
	      vof1=vof1/2; vof2=vof2/1; vof3=vof3/2; vof4=vof4/4;
	    }
	  
	  if ((vof1<=0.5 && vof2>=0.5) || (vof1>=0.5 && vof2<=0.5))
	    {
	      mymatrix[j][i].interplow[0]=x1[0]+(0.5-vof1)*(x2[0]-x1[0])/(vof2-vof1);	      mymatrix[j][i].interplow[1]=x1[1];	      mymatrix[j][i].niptype+=1;	      nip++;
	    }
	  if ((vof4<=0.5 && vof3>=0.5) || (vof4>=0.5 && vof3<=0.5))
	    {
	      if(nip==1)	{
		mymatrix[j][i].interphigh[0]=x4[0]+(0.5-vof4)*(x3[0]-x4[0])/(vof3-vof4);		mymatrix[j][i].interphigh[1]=x4[1];
	      }
	      else 	{
		mymatrix[j][i].interplow[0]=x4[0]+(0.5-vof4)*(x3[0]-x4[0])/(vof3-vof4);		mymatrix[j][i].interplow[1]=x4[1];
	      }
	      nip++;	      mymatrix[j][i].niptype+=5;
	    }
	  if ((vof1<=0.5 && vof4>=0.5) || (vof1>=0.5 && vof4<=0.5))
	    {
	      if (nip==0){
		mymatrix[j][i].interplow[1]=x1[1]+(0.5-vof1)*(x4[1]-x1[1])/(vof4-vof1);		mymatrix[j][i].interplow[0]=x1[0];
	      }
	      else{
		mymatrix[j][i].interphigh[1]=x1[1]+(0.5-vof1)*(x4[1]-x1[1])/(vof4-vof1);		mymatrix[j][i].interphigh[0]=x1[0];
	      }
	      nip++;	      mymatrix[j][i].niptype+=9;
	    }
	  if ((vof3<=0.5 && vof2>=0.5) || (vof3>=0.5 && vof2<=0.5))
	    {
	      if (nip==0){
		mymatrix[j][i].interplow[1]=x3[1]+(0.5-vof3)*(x2[1]-x3[1])/(vof2-vof3);		mymatrix[j][i].interplow[0]=x2[0];
	      }
	      else{
		mymatrix[j][i].interphigh[1]=x3[1]+(0.5-vof3)*(x2[1]-x3[1])/(vof2-vof3);		mymatrix[j][i].interphigh[0]=x2[0];
	      }
	      nip++;	      mymatrix[j][i].niptype+=3;
	    }
	  mymatrix[j][i].lenghtofinterface = sqrt(pow(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0],2) + pow(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1],2));
	  i=i-1;
	  if(i==-1)
	    {
	      i=332;
	      j=j+1;
	    }
	}
      end_c_loop(c,t)
	}
}

DEFINE_EXECUTE_AT_END(IntermediateValueGenerator)/* creation of c hat and gamma hat*/
{
  Domain *d;  Thread *t;  cell_t c;  d=Get_Domain(2); FILE *fp;  fp=fopen("InterValGen.txt","a"); int count=0;
  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
	{
	  if (C_VOF(c,t)!=0 && C_VOF(c,t)!=1)
	    {
	      C_UDSI(c,t,1)=gaminfy*C_UDSI(c,t,0)/(b+C_UDSI(c,t,0));	 
	      fprintf(fp,"bulk conc in this cell is : %g        surface excess in this cell is : %g\n",C_UDSI(c,t,0),C_UDSI(c,t,1));
	    }
	}
      end_c_loop(c,t)
	}
 
  fclose(fp);
}
DEFINE_EXECUTE_AT_END(SourceTerms)/*Calculates the subterms of the source term and stores in the mymatrix*/
{
  Domain *d;  Thread *t,*t1;  cell_t c,c1;  real xc[ND_ND],xc1[ND_ND]; real N1[ND_ND], N2[ND_ND];real mag1,mag2;real ux,uy;
  d=Get_Domain(1);
  int i=332,j=0; 
  real diffusion_term=0;
  real normal_movement=0;
  real compression=0;
  real gammasurfgrad[ND_ND];
  real dilation, marongani;
  FILE *fp;  fp=fopen("SourceTerms.txt","a");
  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
	{
	  C_CENTROID(xc,c,t);
	  diffusion_term=0;normal_movement=0;compression=0;
	  gammasurfgrad[0]=0;
	  gammasurfgrad[1]=0;
	  dilation=0;
	  marongani=0;
	  mymatrix[j][i].source=0;
	  if (mymatrix[j][i].niptype==4)
	    {
	      mag1=sqrt(pow(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0],2) + pow(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1],2));	     
	      N1[0]=(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0])/mag1;	      N1[1]=(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1])/mag1;
	      mag2=sqrt(pow(-mymatrix[j][i].interplow[0]+mymatrix[j][i].interphigh[0],2) + pow(-mymatrix[j][i].interplow[1]+mymatrix[j][i].interphigh[1],2));	     
	      N2[0] =( -mymatrix[j][i].interplow[0] + mymatrix[j][i].interphigh[0])/mag2;	      N2[1] =( -mymatrix[j][i].interplow[1] + mymatrix[j][i].interphigh[1])/mag2;
	     
	      j=j+1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc1[1]-xc[1]);	      j=j-1;
	      gammasurfgrad[1] = C_UDSI(mymatrix[j+1][i].c,mymatrix[j+1][i].t,1)/(xc1[1]-xc[1])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc1[1]-xc[1]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;
	      compression+=(ux*N1[0] + uy*N1[1])*C_UDSI(c,t,1);
	      
	      i=i+1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc1[0]-xc[0]);	      i=i-1;
	      gammasurfgrad[0]=C_UDSI(mymatrix[j][i+1].c,mymatrix[j][i+1].t,1)/(xc1[0]-xc[0])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc1[0]-xc[0]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;
	      compression+=(ux*N2[0] + uy*N2[1])*C_UDSI(c,t,1);
	     
	      normal_movement+=C_U(c,t)*C_UDSI(c,t,1);
	      normal_movement+=C_V(c,t)*C_UDSI(c,t,1);
	    }
	  if (mymatrix[j][i].niptype==6)
	    {
	      mag1=sqrt(pow(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0],2) + pow(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1],2));	     
	      N1[0]=(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0])/mag1;	      N1[1]=(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1])/mag1;
	      mag2=sqrt(pow(-mymatrix[j][i].interplow[0]+mymatrix[j][i].interphigh[0],2) + pow(-mymatrix[j][i].interplow[1]+mymatrix[j][i].interphigh[1],2));	     
	      N2[0] =( -mymatrix[j][i].interplow[0] + mymatrix[j][i].interphigh[0])/mag2;	      N2[1] =( -mymatrix[j][i].interplow[1] + mymatrix[j][i].interphigh[1])/mag2;

	      j=j+1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc1[1]-xc[1]);	      j=j-1;
	      gammasurfgrad[1] = C_UDSI(mymatrix[j+1][i].c,mymatrix[j+1][i].t,1)/(xc1[1]-xc[1])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc1[1]-xc[1]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;
	      compression+=(ux*N1[0] + uy*N1[1])*C_UDSI(c,t,1);
	      
	      j=j-1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc[1]-xc1[1]);	      j=j+1;
	      gammasurfgrad[1]  = gammasurfgrad[1] + C_UDSI(mymatrix[j-1][i].c,mymatrix[j-1][i].t,1)/(xc[1]-xc1[1])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc[1]-xc1[1]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;
	      compression+=(ux*N2[0] + uy*N2[1])*C_UDSI(c,t,1);
	    }
	  if (mymatrix[j][i].niptype==10)
	    {
	      mag1=sqrt(pow(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0],2) + pow(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1],2));	     
	      N1[0]=(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0])/mag1;	      N1[1]=(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1])/mag1;
	      mag2=sqrt(pow(-mymatrix[j][i].interplow[0]+mymatrix[j][i].interphigh[0],2) + pow(-mymatrix[j][i].interplow[1]+mymatrix[j][i].interphigh[1],2));	     
	      N2[0] =( -mymatrix[j][i].interplow[0] + mymatrix[j][i].interphigh[0])/mag2;	      N2[1] =( -mymatrix[j][i].interplow[1] + mymatrix[j][i].interphigh[1])/mag2;     
	     
	      j=j+1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc1[1]-xc[1]);	      j=j-1;
	       gammasurfgrad[1] = C_UDSI(mymatrix[j+1][i].c,mymatrix[j+1][i].t,1)/(xc1[1]-xc[1])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc1[1]-xc[1]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;	      compression+=(ux*N1[0] + uy*N1[1])*C_UDSI(c,t,1);
	     
	      i=i-1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc1[0]-xc[0]);	      i=i+1;
	      gammasurfgrad[0] = C_UDSI(mymatrix[j][i-1].c,mymatrix[j][i-1].t,1)/(xc[0]-xc1[0])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc[0]-xc1[0]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;	      compression+=(ux*N2[0] + uy*N2[1])*C_UDSI(c,t,1);
	     
	      normal_movement+=-C_U(c,t)*C_UDSI(c,t,1);	      normal_movement+=C_V(c,t)*C_UDSI(c,t,1);
	    }
	  if (mymatrix[j][i].niptype==8)
	    {
	      mag1=sqrt(pow(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0],2) + pow(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1],2));	     
	      N1[0]=(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0])/mag1;	      N1[1]=(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1])/mag1;
	      mag2=sqrt(pow(-mymatrix[j][i].interplow[0]+mymatrix[j][i].interphigh[0],2) + pow(-mymatrix[j][i].interplow[1]+mymatrix[j][i].interphigh[1],2));	     
	      N2[0] =( -mymatrix[j][i].interplow[0] + mymatrix[j][i].interphigh[0])/mag2;	      N2[1] =( -mymatrix[j][i].interplow[1] + mymatrix[j][i].interphigh[1])/mag2;	     	     
	     
	      j=j-1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc[1]-xc1[1]);	      j=j+1;
	      gammasurfgrad[1] = C_UDSI(mymatrix[j-1][i].c,mymatrix[j-1][i].t,1)/(xc[1]-xc1[1])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc[1]-xc1[1]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;	      compression+=(ux*N1[0] + uy*N1[1])*C_UDSI(c,t,1);
	     
	      i=i+1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc1[0]-xc[0]);	      i=i-1;
	       gammasurfgrad[0] = C_UDSI(mymatrix[j][i+1].c,mymatrix[j][i+1].t,1)/(xc1[0]-xc[0])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc1[0]-xc[0]);
	     
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;	      compression+=(ux*N2[0] + uy*N2[1])*C_UDSI(c,t,1);
	     
	      normal_movement+=C_U(c,t)*C_UDSI(c,t,1);	      normal_movement+=-C_V(c,t)*C_UDSI(c,t,1);
	    }
	  if (mymatrix[j][i].niptype==12)
	    {
	      mag1=sqrt(pow(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0],2) + pow(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1],2));	     
	      N1[0]=(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0])/mag1;	      N1[1]=(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1])/mag1;
	      mag2=sqrt(pow(-mymatrix[j][i].interplow[0]+mymatrix[j][i].interphigh[0],2) + pow(-mymatrix[j][i].interplow[1]+mymatrix[j][i].interphigh[1],2));	     
	      N2[0] =( -mymatrix[j][i].interplow[0] + mymatrix[j][i].interphigh[0])/mag2;	      N2[1] =( -mymatrix[j][i].interplow[1] + mymatrix[j][i].interphigh[1])/mag2;	     	     
	     
	      i=i+1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc1[0]-xc[0]);	      i=i-1;
	      gammasurfgrad[0] = C_UDSI(mymatrix[j][i+1].c,mymatrix[j][i+1].t,1)/(xc1[0]-xc[0]) - C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc1[0]-xc[0]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;	      compression+=(ux*N1[0] + uy*N1[1])*C_UDSI(c,t,1);
	     
	      i=i-1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc[0]-xc1[0]);	      i=i+1;
	      gammasurfgrad[0] = gammasurfgrad[0] - C_UDSI(mymatrix[j][i-1].c,mymatrix[j][i-1].t,1)/(xc[0]-xc1[0]) - C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc[0]-xc1[0]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;	      compression+=(ux*N2[0] + uy*N2[1])*C_UDSI(c,t,1);
	    }
	  if (mymatrix[j][i].niptype==14)
	    {
	      mag1=sqrt(pow(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0],2) + pow(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1],2));	     
	      N1[0]=(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0])/mag1;	      N1[1]=(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1])/mag1;
	      mag2=sqrt(pow(-mymatrix[j][i].interplow[0]+mymatrix[j][i].interphigh[0],2) + pow(-mymatrix[j][i].interplow[1]+mymatrix[j][i].interphigh[1],2));	     
	      N2[0] =( -mymatrix[j][i].interplow[0] + mymatrix[j][i].interphigh[0])/mag2;	      N2[1] =( -mymatrix[j][i].interplow[1] + mymatrix[j][i].interphigh[1])/mag2;	     	     
	     
	      j=j-1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc[1]-xc1[1]);	      j=j+1;
	      gammasurfgrad[1] = C_UDSI(mymatrix[j-1][i].c,mymatrix[j-1][i].t,1)/(xc[1]-xc1[1])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc[1]-xc1[1]);
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;	      compression+=(ux*N1[0] + uy*N1[1])*C_UDSI(c,t,1);
	     
	      i=i-1;	      c1=mymatrix[j][i].c;	      t1=mymatrix[j][i].t;	      C_CENTROID(xc1,c1,t1);	      diffusion_term=diffusion_term+(D-Di)*(C_UDSI(c1,t1,1)-C_UDSI(c,t,1))/(xc[0]-xc1[0]);	      i=i+1;
	      gammasurfgrad[0] = C_UDSI(mymatrix[j][i-1].c,mymatrix[j][i-1].t,1)/(xc[0]-xc1[0])-C_UDSI(mymatrix[j][i].c,mymatrix[j][i].t,1)/(xc[0]-xc1[0]);	  
	      ux=C_U(c1,t1)/2+C_U(c,t)/2;	      uy=C_V(c1,t1)/2+C_V(c,t)/2;	      compression+=(ux*N2[0] + uy*N2[1])*C_UDSI(c,t,1);
	     
	      normal_movement+=-C_U(c,t)*C_UDSI(c,t,1);	      normal_movement+=-C_V(c,t)*C_UDSI(c,t,1);
	    }
	  dilation = (C_U(c,t)*gammasurfgrad[0] + C_V(c,t)*gammasurfgrad[1])*mymatrix[j][i].lenghtofinterface;
	  marongani = 2*(b + C_UDSI(c,t,0))*D*(pow(gammasurfgrad[0],2) + pow(gammasurfgrad[1],2))*mymatrix[j][i].lenghtofinterface/(gaminfy*b);
	  mymatrix[j][i].source = diffusion_term + normal_movement - compression - marongani + dilation;
	  /*	  fprintf(fp,"Cell centroid : (%g  ,  %g)\n",xc[0],xc[1]);
		  fprintf(fp," diffusion_term: %g          normal_movement: %g         compression : %g       surface gradient of gamma :%g      dilation = %g     marongani = %g\n\n",-diffusion_term,normal_movement,-compression,sqrt(pow(gammasurfgrad[0],2) + pow(gammasurfgrad[1],2)),dilation, -marongani);*/
	  i=i-1;
	  if(i==-1)
	    {
	      i=332;
	      j=j+1;
	    }	   
	}
      end_c_loop(c,t)
	}
  fclose(fp);
}

DEFINE_EXECUTE_AT_END(finalUpdater)
{
  Domain *d;  Thread *t;  cell_t c;  real xc[ND_ND];  d=Get_Domain(2);  real conc, gamma,conchat,gammahat;  real b=0.000663;  real gaminfy=0.00000291;  real lenghtofinterface, Nhat, N, source;
  real A,B,C;
  int i=332,j=0;  FILE *fp; fp=fopen("finalupdater.txt","a");
  real deltaT;
  deltaT = CURRENT_TIMESTEP;
  printf("time step is : %g", deltaT);
  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
	{
	  C_CENTROID(xc,c,t);
	  if (C_VOF(c,t)!=0 && C_VOF(c,t)!=1)
	    {
	      conchat=C_UDSI(c,t,0);
	      gammahat=C_UDSI(c,t,1);  fprintf(fp,"VOF : %g    bulk conc: %g     surface excess: %g     interval gen values",C_VOF(c,t),conchat,gammahat);
	      lenghtofinterface=sqrt(pow(mymatrix[j][i].interplow[0]-mymatrix[j][i].interphigh[0],2) + pow(mymatrix[j][i].interplow[1]-mymatrix[j][i].interphigh[1],2));
	      fprintf(fp,"Lenght of interface : %g   ",lenghtofinterface);
	      Nhat=conchat*C_VOLUME(c,t)*C_VOF(c,t) + gammahat*lenghtofinterface;	   fprintf(fp,"Nhat : %g    ",Nhat);    
	      source=mymatrix[j][i].source;
	      N=source*deltaT + Nhat;	     fprintf(fp,"N : %g    ",N);
	      
	      B=b*C_VOLUME(c,t)*C_VOF(c,t) + gaminfy*lenghtofinterface - N;
	      A=C_VOF(c,t)*C_VOLUME(c,t);
	      C=b*N;
 	      conc=(-B+sqrt(B*B+4*C*A))/(2*A);	      C_UDSI(c,t,0)=conc;	      gamma=gaminfy*conc/(b+conc);	      C_UDSI(c,t,1)=gamma;
	      fprintf(fp,"bulk conc in this cell is : %g          surface excess in this cell : %g\n",C_UDSI(c,t,0),C_UDSI(c,t,1));
	      /* fprintf(fp,"Centroid : (%g  ,  %g),      Lenght of interface : %g              source : %g      Nhat: %g \n",xc[0],xc[1],lenghtofinterface,source,Nhat);
              if(lenghtofinterface==0){fprintf(fp,"Interplow : (%g,     %g),        Interphigh : (  %g,     %g)\n\n", mymatrix[j][i].interplow[0],mymatrix[j][i].interplow[1], mymatrix[j][i].interphigh[0],mymatrix[j][i].interphigh[1]);}
	      */
	      }
	  i=i-1;
	  if(i==-1)
	    {
	      i=332;
	      j=j+1;
	    }
	}
      end_c_loop(c,t)
	}
  fclose(fp);
}
DEFINE_EXECUTE_AT_END(massBalanceCheck)
{
  Domain *d;  Thread *t;  float gammatotal=0,gammaair=0,gammawater=0,gammainter=0, conctotal=0,concwater=0,concair=0,concinter=0;    cell_t c;   d=Get_Domain(2);
  FILE *fp;
  fp=fopen("Mass at all times.txt","a");
  int i=332,j=0; 
  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
	{
	  gammatotal = gammatotal + mymatrix[j][i].lenghtofinterface * C_UDSI(c,t,1);
	  conctotal = conctotal + C_VOLUME(c,t) * C_UDSI(c,t,0) * C_VOF(c,t);
	  if(C_VOF(c,t)==0){ gammaair = gammaair + mymatrix[j][i].lenghtofinterface * C_UDSI(c,t,1);	  concair = concair + C_VOLUME(c,t) * C_UDSI(c,t,0) * C_VOF(c,t);}
	  if(C_VOF(c,t)==1){ gammawater = gammawater + mymatrix[j][i].lenghtofinterface * C_UDSI(c,t,1);	  concwater = concwater + C_VOLUME(c,t) * C_UDSI(c,t,0) * C_VOF(c,t);}
	  if(C_VOF(c,t)!=1 && C_VOF(c,t)!=0){ gammainter = gammainter + mymatrix[j][i].lenghtofinterface * C_UDSI(c,t,1);	  concinter = concinter + C_VOLUME(c,t) * C_UDSI(c,t,0) * C_VOF(c,t);}
	}
      end_c_loop(c,t)
	}
  i=i-1;
  if(i==-1)
    {
      i=332;
      j=j+1;
    }
  real time;
  time=CURRENT_TIME;
  fprintf(fp,"gammainter : %g, gammaair : %g,  gammawater : %g,  gammatotal: %g\nconcinter : %g,   concwater : %g,  concair : %g,   conctotal : %g,    time: %g\n\n", gammainter,gammaair,gammawater,gammatotal,concinter,concwater,concair,conctotal, time);
  fclose(fp);
}
DEFINE_EXECUTE_AT_END(velTAverageAndGammaDistribution)
{
  Domain *d;  Thread *t;  float sum_velx=0;  float count=0;  cell_t c;  real xc[ND_ND];  d=Get_Domain(2);
  real gammatotal=0,gammawater=0,gammaair=0,gammainter=0,conctotal=0,concair=0,concwater=0,concinter=0;
  FILE *fp;
  fp=fopen("Gamma distribution and velocity.txt","a");
  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
	{
	  C_CENTROID(xc,c,t);
	  gammatotal = gammatotal + C_UDSI(c,t,1);
	  conctotal = conctotal + C_UDSI(c,t,0);
	  if(C_VOF(c,t)==0){ gammaair = gammaair + C_UDSI(c,t,1);	  concair = concair + C_UDSI(c,t,0);}
	  if(C_VOF(c,t)==1){ gammawater = gammawater + C_UDSI(c,t,1);	  concwater = concwater + C_UDSI(c,t,0);}
	  if(C_VOF(c,t)!=1 && C_VOF(c,t)!=0){ gammainter = gammainter + C_UDSI(c,t,1);	  concinter = concinter + C_UDSI(c,t,0);}
	  if(C_VOF(c,t)!=0 && C_VOF(c,t)!=1)
	    {
	      count=count+1;
	      sum_velx=sum_velx+C_U(c,t);
	    }
	}
      end_c_loop(c,t)
	}
  real avgx=sum_velx/count;
  real time;
  time=CURRENT_TIME;
  fprintf(fp,"gammainter : %g, gammaair : %g,  gammawater : %g,  gammatotal: %g\nconcinter : %g,   concwater : %g,  concair : %g,   conctotal : %g\n", gammainter,gammaair,gammawater,gammatotal,concinter,concwater,concair,conctotal);
  fprintf(fp,"bubbleVelocity: %g, time: %g\n\n",avgx,time);
  fclose(fp);
}

DEFINE_EXECUTE_AT_END(xaspectratiocalculator)
{
  Domain *d;  Thread *t;  cell_t c;  real xc[ND_ND];  d=Get_Domain(2);
  real x_low=0.11,x_high=0,y_low=0.11,y_high=0;
  real aspectratio;
  FILE *fp;
  fp=fopen("aspect ratio.txt","a");
  thread_loop_c(t,d)
    {
      begin_c_loop(c,t)
	{
	  if(C_VOF(c,t)!=0 && C_VOF(c,t)!=1)
	    {
	      C_CENTROID(xc,c,t);
	      if(x_low>xc[0]){x_low=xc[0];}
	      if(x_high<xc[0]){x_high=xc[0];}
	      if(y_low>xc[1]){y_low=xc[1];}
	      if(y_high<xc[1]){y_high=xc[1];}
	    }
	}
      end_c_loop(c,t)
	}
  aspectratio = (y_high-y_low)/(x_high-x_low);
  real time;
  time=CURRENT_TIME;
  fprintf(fp,"aspect ratio : %g,   time : %g\n",aspectratio,time);
  fclose(fp);
}

DEFINE_PROPERTY(SurfaceTension,c,t)
{
  real surfaceTension, a=0.01,gaminfy= 0.00000291;
  if(C_UDSI_M1(c,t,1)<=gaminfy)
    {
      surfaceTension=0.03-a*C_UDSI_M1(c,t,1)/gaminfy;
    }
  else 
    {
      surfaceTension=0.03-a;
    }
  return surfaceTension;
}


