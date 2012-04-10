/*
 * =====================================================================================
 *
 *       Filename:  ux_point.c
 *
 *    Description: use discontinuous Galerkin,polynomial truncation for the nonlinear  terms in space,and Runge-kutta in time,to compute the conserfvation law
 *    equation u_t+u_x=0 with inital codition u(x,0)=sin(\pi,x) and periodic 
 *    boundary condition
 *
 *
 *        Version:  1.0
 *        Created:  2011Äê11ÔÂ05ÈÕ 10Ê±13·Ö49Ãë
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define	max(A,B)	((A)>(B) ? (A) : (B))
#define min(A,B)	((A)<(B) ? (A) : (B))
#define order 5
static double ma[5];
//double error[3];
double dtime;
double fle(int mp,double x);
double fled(int mp,double x);
//double min(double* x1,double* x2);
//double max(double* x1,double* x2);
double eval(double* a,double x);
int mp;
int	mo;
double flux(double x);
void res(double **u,double **hg, double *x,double *dx,int n);
  void rk(double **u,double **hg,double *x,double *dx,int n);
void  compute(double **u,double **hg,double *x,double *dx,int n);


int main(int argc, char* argv[])
{
//	double er1,er2,er3;
//	double rr1,rr2,rr3;
	double **u;
	double **hg;
	double *x;
	double *dx;
	//L_1 error
	double error1[6];
	//L_2 error
	double error2[6];
	//L_{infty} error
	char name[12];
	double  error3[6];
	double temp1,temp2,temp3;
	double *test;
	int i,j,k;
	int n=10;
	FILE *fp;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	//ÖÊÁ¿ŸØÕóµÄÄæŸØÕó
	ma[0]=1.0;
	ma[1]=12.0;
	ma[2]=180.0;
	ma[3]=2800.0;
	ma[4]=44100.0;
	printf("input the order(mo=1,2,3,4,5):\n");
	scanf("%d",&mo);
	mp=mo-1;
	//ÇøŒä·Ö±ðÆÊ·Ö10,20,40,80,160·Ý
	for(j=0;j<6;j++)
	{
		//·ÖÅäÄÚŽæ
		//ÏÈ·ÖÅäÖžÕëµ¥Ôª
		u=(double **)malloc(n*sizeof(double *));
		hg=(double **)malloc(n*sizeof(double *));
		//ÔÚ·ÖÅä×Ö·ûµ¥Ôª
		for(i=0;i<n;i++)
		{
			u[i]=(double *)malloc(mo*sizeof(double));
			hg[i]=(double *)malloc(mo*sizeof(double));
		}

		x=(double *)malloc(n*sizeof(double));
		dx=(double *)malloc(n*sizeof(double));		

		compute(u,hg,x,dx,n);
		printf("%d intervals has complete, now we tackle %d intevals \n",n,n*2);
		//ÊÍ·ÅÄÚŽæ 
		for(i=0;i<n;i++)
		{
			free((void *)u[i]);
			free((void *)hg[i]);
		}
		free((void *)u);
		free((void *)hg);
		free(x);
		free(dx);
		n=n*2;
		//printf("%d\n\n",n);
	}


	//ŒÆËãÎó²îœ×
	//¶ÁÈëÎó²îÊýŸÝ
	sprintf(name,"%d%s",mo,"error");
	fp=fopen(name,"r+");
	i=0;
	while((read=getline(&line,&len,fp))!=-1)
	{
	sscanf(line,"%lf\t%lf\t%lf",&error1[i],&error2[i],&error3[i]);
	i++;
	}
	//¿ªÊŒŒÆËãÎó²î
	k=i-1;
	printf("Error:\n");
	printf("L_1\tL_2\tL_{infty}\n");
	for(i=0;i<k;i++)
	{
	temp1=log(error1[i]/error1[i+1])/log(2.0);
	temp2=log(error2[i]/error2[i+1])/log(2.0);
	temp3=log(error3[i]/error3[i+1])/log(2.0);
	printf("%5.4f\t%5.4f\t%5.4f\n",temp1,temp2,temp3);
	}
	return 0; 
}

void  compute(double **u,double **hg, double *x,double *dx,int n){
	//the order in space of the scheme 
	//the cflc number
	double cflc;
	/*the terminal mytime:tprint;
	*the tprint can't be too large since the exact solution is e^{-t \pi^2}
	* *sin(x\pi),If t is too large ,the solution will be very very small.
	*/
	// the maximum number of mytime steps:kcmax
	//int kcmax = 1000000;
	double tprint=2.0;
	double xleft=0.0,xright=2.0;
	double xlen,dxuni;
	double dxmin;
	//double temp[n+1];
	double *temp;
	int i,j;
	//double a[mo];
	double mytime;
	//double bl[mo],br[mo];
	double pi=4.0*atan(1.0);
	double yy;
	double xvalue;
	double uexact;
	double unum;
	double rr;
	double error1=0.0;
	double error2=0.0;
	double error3=0.0;
	//int t=0;
	FILE* fp;
	FILE* ep;
	char name[12];
	char ename[12];
	//intial data
	xlen=xright-xleft;
	dxuni=xlen/n;
	dxmin = dxuni;
	temp=(double *)malloc((n+1)*sizeof(double));
	for(i = 0;i<n+1;i++)
		temp[i]=xleft+dxuni*i;
	//ÉèÖÃËæ»úÊýÉú³ÉÖÖ×Ó
	//srand(time(NULL));
	//¶ÔÇøŒäœøÐÐ²»ŸùÔÈÆÊ·Ý
//	for(i=1;i<n;i++)
//		temp[i]=temp[i]+0.1*dxuni*(((double)rand())/RAND_MAX-0.5);
	//¶ÔÅŒœÚµã
	for(i=0;i<n;i++)
	{
		x[i]=(temp[i+1]+temp[i])/2.0;
		dx[i]=temp[i+1]-temp[i];
	}
	for(i=0;i<n;i++)
		for(j=0;j<mo;j++)
			u[i][j]=ma[j]*((sin(pi*(x[i]+(-0.5)*dx[i]))*fle(j,-0.5)+sin(pi*(x[i]+0.5*dx[i]))*fle(j,0.5))*(1.0/30.0)+(sin(pi*(x[i]+(-0.142615758)*dx[i]))*fle(j,-0.142615758)+sin(pi*(x[i]+0.142615758*dx[i]))*fle(j,0.142615758))*0.277429189+(sin(pi*(x[i]+(-0.382527662)*dx[i]))*fle(j,-0.382527662)+sin(pi*(x[i]+0.382527662*dx[i]))*fle(j,0.382527662))*0.189237478);

	sprintf(name,"%d%d",mo,n);
	fp=fopen(name,"a+");
	fprintf(fp,"%d intervals 		order:%d\n",n,mo);

	//determain the cflc number
	if (mp==1)
		cflc = 0.3;
	else if (mp==2)
		cflc = 0.18;
	else if (mp==3)
		cflc = 0.1;
	else 
		cflc = 0.08;
	mytime=0.0;
	if(mp < 2 + 1.0e-5)
		rr=1.0;
	else
		rr=(mp+1.0)/4.0;
	dtime=cflc*pow(dxmin,rr);
	//printf("%f\n",dtime);
	if(mytime+dtime > tprint)
		dtime=tprint-mytime;
	while(mytime <= tprint)
	{
		rk(u,hg,x,dx,n);
		mytime+=dtime;
	}
	
//	mytime = mytime-dtime;
	fprintf(fp,"coefficient£º\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<mo;j++)
		{
			fprintf(fp,"%f  ",u[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"x's value\tnumerical solution\texact solution\n");

	//compute the error error[0]=L_1,error[1]=L_2,error[3]=L_{\infty}
/* for(i=0;i<3;i++)
 * 		error[i]=0.0;
 */
	for(i=0;i<n;i++)
		for(j=0;j<41;j++)
		{
			if(j==20)
			{
				xvalue=x[i];
				unum=eval(u[i],0);
				uexact=sin(pi*(x[i]-mytime));
				fprintf(fp,"%12.8f\t%12.8f\t%12.8f\n",xvalue,unum,uexact);
				yy=fabs(unum-uexact);
			}
			else
			{
				yy=fabs(eval(u[i],(j-20.0)/40.0)-sin(pi*(x[i]+((j-20.0)/40.0)*dx[i]-mytime)));
			}
			error1=error1+yy*dx[i];
			error2=error2+yy*yy*dx[i];
			error3=max(error3,yy);

		}
	sprintf(ename,"%d%s",mo,"error");
	ep=fopen(ename,"a+");

	error1=error1/41.0/2.0;
	error2=sqrt(error2/41.0/2.0);
	//fprintf(fp,"error:L_1\tL_2\tL_{infty}\n");
	fprintf(fp,"%20.18f\t%20.18f\t%20.18f\n",error1,error2,error3);
	fprintf(ep,"%20.18f\t%20.18f\t%20.18f\n",error1,error2,error3);
	
	free(temp);
	fclose(fp);
	fclose(ep);	
}
//The function of Legendre polynomial
double fle(int k,double x)
{
	
	double value;
	if(k==0)
		value=1.0;
	else if(k==1)
		value = x;
	else if(k==2)
		value = x*x-1.0/12.0;
	else if(k==3)
		value = x*x*x-0.15*x;
	else if(k==4)
		value = (x*x-3.0/14.0)*x*x+3.0/560.0;
	else
		value = 0.0;
	return value;	
}
//the function of derivative of Legendre polynomial
double fled(int k,double x)
{
	double value;
	if (k==0)
		value = 0.0;
	else if(k==1)
		value = 1.0;
	else if(k==2)
		value = 2.0*x;
	else if(k==3)
		value = 3.0*x*x-0.15;
	else if(k==4)
		value = (4.0*x*x-3.0/7.0)*x;
	else
		value = 0.0;
	return value;
}
/*
double min(double *x1,double *x2)
{
	if(*x1 < *x2)
		return *x1;
	else
		return *x2;
}
*/
double eval(double *a,double x)
{
	int i;
	double value=0.0;
	for(i=0;i<mo;i++)
	{
		value += a[i]*fle(i,x);
	}
		return value;
}
// the function of flux
double flux(double x)
{
	return x;
}
void rk(double **u,double **hg,double *x,double *dx,int n)
{
	int i,j;
	double **v;
	double **temp1;
	double **temp2;
	double **temp3;
	/*
	double v[n][mo];
	double temp1[n][mo];
	double temp2[n][mo];
	double temp3[n][mo];
	*/

	//·ÖÅäÄÚŽæ
		//ÏÈ·ÖÅäÖžÕëµ¥Ôª
		v=(double **)malloc(n*sizeof(double *));
		temp1=(double **)malloc(n*sizeof(double *));
		temp2=(double **)malloc(n*sizeof(double *));
		temp3=(double **)malloc(n*sizeof(double *));
		//ÔÚ·ÖÅä×Ö·ûµ¥Ôª
		for(i=0;i<n;i++)
		{
			v[i]=(double *)malloc(mo*sizeof(double));
			temp1[i]=(double *)malloc(mo*sizeof(double));
			temp2[i]=(double *)malloc(mo*sizeof(double));
			temp3[i]=(double *)malloc(mo*sizeof(double));
		}


	
	for(i=0;i<n;i++)
		for(j=0; j<mo; j++){
			v[i][j]=u[i][j];
		}
	res(u,hg,x,dx,n);
	for(i=0;i<n;i++)
		for(j=0;j<mo;j++)
		{
			u[i][j]=v[i][j]+0.5*dtime*hg[i][j];
			temp1[i][j]=u[i][j];
		}

	res(u,hg,x,dx,n);
	for(i=0;i<n;i++)
		for(j=0;j<mo;j++)
		{
			u[i][j]=v[i][j]+0.5*dtime*hg[i][j];
			temp2[i][j]=u[i][j];
		}

	res(u,hg,x,dx,n);
	for(i=0;i<n;i++)
		for(j=0;j<mo;j++)
		{
			u[i][j]=v[i][j]+dtime*hg[i][j];
			temp3[i][j]=u[i][j];
		}

	res(u,hg,x,dx,n);
	for(i=0;i<n;i++)
		for(j=0;j<mo;j++)
			u[i][j]=(-v[i][j]+temp1[i][j]+2*temp2[i][j]+temp3[i][j]+
					0.5*dtime*hg[i][j])/3.0;
	//ÊÍ·ÅÄÚŽæ 
		for(i=0;i<n;i++)
		{
			free((void *)v[i]);
			free((void *)temp1[i]);
			free((void *)temp2[i]);
			free((void *)temp3[i]);
		}
		free((void *)v);
		free((void *)temp1);
		free((void *)temp2);
		free((void *)temp3);

}
void res(double **u,double **hg,double *x,double *dx,int n)
{
	/*
	double h[n][mo];
	double un[n];
	double unn[n];
	double up[n];
	double bl[mo];
	double br[mo];
	*/
	int i,j,k;
	double **h;
	double *un;
	double *unn;
	double *up;
	double *bl;
	double *br;
	//·ÖÅäÄÚŽæ
		//ÏÈ·ÖÅäÖžÕëµ¥Ôª
		h=(double **)malloc(n*sizeof(double *));

		//ÔÚ·ÖÅä×Ö·ûµ¥Ôª
		for(i=0;i<n;i++)
		{
			h[i]=(double *)malloc(mo*sizeof(double));			
		}
		un=(double *)malloc(n*sizeof(double));
		unn=(double *)malloc(n*sizeof(double));
		up=(double *)malloc(n*sizeof(double));
		bl=(double *)malloc(mo*sizeof(double));
		br=(double *)malloc(mo*sizeof(double));
	for(i=0;i<mo;i++)
	{
		bl[i]=fle(i,-0.5);
		br[i]=fle(i,0.5);
	}
	for(i=0;i<n;i++)
	{
		un[i]=eval(u[i],0.5);
		up[i]=eval(u[i],-0.5);
	}
	unn[0]=un[n-1];
	for(i=1;i<n;i++)
		unn[i]=un[i-1];
	for(i=0;i<n;i++)
		for(j=0;j<mo;j++)
			h[i][j]=-un[i]*br[j]+unn[i]*bl[j];
	if(mp > 0)
			for(i=0;i<n;i++)
			{
				if(mp==1)
					h[i][1]=h[i][1]+(flux(eval(u[i],-0.5))+
							flux(eval(u[i],0.0))*4.0+flux(eval(u[i],0.5)))/6.0;
//				else if(mp==2)
//					for(k=1;k<mo;k++)
//						h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*(1.0/12.0)+
//							(flux(eval(u[i],-sqrt(5.0)/10.0))*fled(k,-sqrt(5.0)/10.0)+
//						flux(eval(u[i],sqrt(5.0)/10.0))*fled(k,sqrt(5.0)/10.0))*(5.0/12.0));
//				else if(mp==3)
//					for(k=1;k<mo;k++)
//						h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*0.05+
//						(flux(eval(u[i],-sqrt(21.0)/14.0))*fled(k,-sqrt(21.0)/14.0)+
//							flux(eval(u[i],sqrt(21.0)/14.0))*fled(k,sqrt(21.0)/14.0))*(49.0/180.0)+
//						(flux(eval(u[i],0))*fle(k,0))*(16.0/45.0));
//
				else
					for(k=1;k<mo;k++)
						h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*(1.0/30.0)+(flux(eval(u[i],sqrt(147.0-42.0*sqrt(7.0))/42.0))*fled(k,sqrt(147.0-42.0*sqrt(7.0))/42.0)+flux(eval(u[i],-sqrt(147.0-42.0*sqrt(7.0))/42.0))*fled(k,-sqrt(147.0-42.0*sqrt(7.0))/42.0))*((14.0+sqrt(7.0))/60.0)+(flux(eval(u[i],sqrt(147.0+42.0*sqrt(7.0))/42.0))*fled(k,sqrt(147.0+42.0*sqrt(7.0))/42.0)+flux(eval(u[i],-sqrt(147.0+42.0*sqrt(7.0))/42.0))*fled(k,-sqrt(147.0+42.0*sqrt(7.0))/42.0))*((14.0-sqrt(7.0))/60.0));
			}
//	for(i=0;i<n;i++)
//		for(k=0;k<mo;k++)
//			hg[i][k]=ma[k]*h[i][k]/dx[i];
////	for(i=0;i<n;i++)
//	{
//		un[i]=eval((*(w+i)),0.5);
//		up[i]=eval((*(w+i)),-0.5);
//	}
//	unn[n-1]=up[0];
//	for(i=0;i<n-1;i++)
//		unn[i]=up[i+1];
//	for(i=0;i<n;i++)
//		for(j=0;j<mo;j++)
//			h[i][j]=unn[i]*br[j]-up[i]*bl[j];
//	if(mp > 0)
//		for(i=0;i<n;i++)
//		{
//			if(mp==1)
//				h[i][1]=h[i][1]+(flux(up[i])+
//						flux(eval((*(w+i)),0.0))*4.0+flux(un[i]))/6.0;
//	/* 		else if(mp==2)
//	 * 			for(k=1;k<mo;k++)
//	 * 				h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*0.083333333+
//	 * 						(eval((*(w+i)),-0.223606798)*fled(k,-0.223606798)+
//	 * 	eval((*(w+i)),0.223606798)*fled(k,0.223606798))*0.416666667);
//	 * 		else if(mp==3)
//	 * 			for(k=1;k<mo;k++)
//	 * 				h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*0.05+(flux(eval((*(w+i)),-0.458257569))*fled(k,-0.458257569)+flux(eval((*(w+i)),0.458257569))*fled(k,0.458257569))*0.272222222+(flux(eval((*(w+i)),0))*fle(k,0))*0.355555556);
//	 */
//			else
//				for(k=1;k<mo;k++)
//					h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*(1.0/30.0)+(flux(eval((*(w+i)),0.142615758))*fled(k,0.142615758)+flux(eval((*(w+i)),-0.142615758))*fled(k,-0.142615758))*0.277429189+(flux(eval((*(w+i)),0.382527662))*fled(k,0.382527662)+flux(eval((*(w+i)),-0.382527662))*fled(k,-0.382527662))*0.189237478);
//		}
	for(i=0;i<n;i++)
		for(k=0;k<mo;k++)
			hg[i][k]=ma[k]*h[i][k]/dx[i];
	//ÊÍ·ÅÄÚŽæ
	for(i=0;i<n;i++)
		free((void *)h[i]);
	free((void *) h);
	free(un);
	free(unn);
	free(up);
	free(bl);
	free(br);
}
/*
double max(double *x1,double *x2)
{
	if(*x1 > *x2)
		return *x1;
	else
		return *x2;
}
*/


