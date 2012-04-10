/*
 * =====================================================================================
 *
 *       Filename:  uxx.c
 *
 *    Description:  
 *
 *
 *        Version:  1.0
 *        Created:  2011年11月05日 10时13分49秒
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
#define order 5
static double ma[5];
//double error[3];
double dtime;
double fle(int mp,double x);
double fled(int mp,double x);
double min(double* x1,double* x2);
double max(double* x1,double* x2);
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
	double error1[5];
	//L_2 error
	double error2[5];
	//L_{infty} error
	double  error3[5];
	double temp1,temp2,temp3;
	int i,j,k;
	int n=10;
	FILE *fp;
	char name[12];
	char *line=NULL;
	size_t len=0;
	ssize_t read;
	//质量矩阵的逆矩阵
	ma[0]=1.0;
	ma[1]=12.0;
	ma[2]=180.0;
	ma[3]=2800.0;
	ma[4]=44100.0;
	printf("input order(mo=1,2,3,4,5):\n");
	scanf("%d",&mo);
	mp=mo-1;
	//区间分别剖分10,20,40,80,160份
	for(j=0;j<4;j++)
	{
		//分配内存
		//先分配指针单元
		u=(double **)malloc(n*sizeof(double *));
		hg=(double **)malloc(n*sizeof(double *));
		//在分配字符单元
		for(i=0;i<n;i++)
		{
			u[i]=(double *)malloc(mo*sizeof(double));
			hg[i]=(double *)malloc(mo*sizeof(double));
		}

		x=malloc(n*sizeof(double));
		dx=malloc(n*sizeof(double));

		compute(u,hg,x,dx,n);
		printf("%d intervals has complete,now we tackle %d intervals \n",n,n*2);
		sleep(2);
		//释放内存 
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
	//计算误差阶
	//读入误差数据
	sprintf(name,"%d%s",mo,"error");
	fp=fopen(name,"r+");
	i=0;
	while((read=getline(&line,&len,fp))!=-1)
	{
		sscanf(line,"%lf\t%lf\t%lf",&error1[i],&error2[i],&error3[i]);
		i++;
	}
	//开始计算误差
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
	double tprint=0.6;
	double xleft=0.0,xright=2.0;
	double xlen,dxuni;
	double dxmin;
	double temp[n+1];
	int i,j;
	//double a[mo];
	double mytime;
	//double bl[mo],br[mo];
	double pi=4.0*atan(1.0);
	double yy;
	double xvalue;
	double uexact;
	double unum;
	double seed;
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
	for(i = 0;i<n+1;i++)
		temp[i]=xleft+dxuni*i;
	//设置随机数生成种子
//	srand(time(NULL));
	//对区间进行不均匀剖份
//	for(i=1;i<n;i++)
//		temp[i]=temp[i]+0.1*dxuni*(((double)rand())/RAND_MAX-0.5);
	//对偶节点
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
	fprintf(fp,"intervals :%d		order:%d\n",n,mo);

	//determain the cflc number
	if (mp==1)
		cflc = 0.06;
	else if (mp==2)
		cflc = 0.002;
	else if (mp==3)
		cflc = 0.004;
	else 
		cflc = 0.002;
	mytime=0.0;
	dtime=cflc*dxmin*dxmin;
	//printf("%f\n",dtime);
	if(mytime+dtime > tprint)
		dtime=tprint-mytime;
	while(mytime <= tprint)
	{
/* 		if(t==1||t==0)
 * 		{
 * 			for(i=0;i<n;i++)
 * 			{
 * 				for(j=0;j<mo;j++)
 * 				{
 * 					fprintf(fp,"%f  ",u[i][j]);
 * 				}
 * 				fprintf(fp,"\n");
 * 			}
 * 		fprintf(fp,"===================================================================================\n");
 * 
 * 
 * 
 * 		}
 */
		rk(u,hg,x,dx,n);
		mytime+=dtime;
	//	printf("%f\n",mytime);
	//	t++;
	}
	
	//mytime = mytime-dtime;
	fprintf(fp,"coeffs：\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<mo;j++)
		{
			fprintf(fp,"%f  ",u[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"x\tnumerical solution\texect solution\n");

	//compute the error error[0]=L_1,error[1]=L_2,error[3]=L_{\infty}
/* for(i=0;i<3;i++)
 * 		error[i]=0.0;
 */
	for(i=0;i<n;i++)
/* 	{
 * 		error[0]=error[0]+((fabs(eval(*(u+i),-0.5)-exp(-mytime*pi*pi)*sin(pi*(x[i]-0.5*dx[i])))+fabs((eval(*(u+i),0.5)-exp(-mytime*pi*pi)*sin(pi*(x[i]+0.5*dx[i])))))*(1.0/30.0)+(fabs(eval(*(u+i),-0.142615758)-exp(-mytime*pi*pi)*sin(pi*(x[i]-0.142615758*dx[i])))+fabs(eval(*(u+i),0.142615758)-exp(-mytime*pi*pi)*sin(pi*(x[i]+0.142615758*dx[i]))))*0.277429189+(fabs(eval(*(u+i),-0.382527662)-exp(-mytime*pi*pi)*sin(pi*(x[i]-0.382527662)*dx[i]))+fabs(eval(*(u+i),0.382527662)-exp(-mytime*pi*pi)*sin(pi*(x[i]+0.382527662)*dx[i])))*0.189237478);
 * 
 * 	}
 */


		for(j=0;j<41;j++)
		{
			if(j==20)
			{
				xvalue=x[i];
				unum=eval((*(u+i)),0);
				uexact=exp(-mytime*pi*pi)*sin(pi*x[i]);
				fprintf(fp,"%12.8f\t%12.8f\t%12.8f\n",xvalue,unum,uexact);
				yy=fabs(unum-uexact);
			}
			else
			{
				yy=fabs(eval((*(u+i)),(j-20.0)/40.0)-exp(-mytime*pi*pi)*sin(pi*(x[i]+((j-20.0)/40.0)*dx[i])));
			}
			error1=error1+yy*dx[i];
			error2=error2+yy*yy*dx[i];
			error3=max(&error3,&yy);

		}
	sprintf(ename,"%d%s",mo,"error");
	ep=fopen(ename,"a+");

	error1=error1/41.0/2.0;
	error2=sqrt(error2/41.0/2.0);
	fprintf(fp,"误差:L_1\tL_2\tL_{infty}\n");
	fprintf(fp,"%20.18f\t%20.18f\t%20.18f\n",error1,error2,error3);
	fprintf(ep,"%20.18f\t%20.18f\t%20.18f\n",error1,error2,error3);

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
		value = (x*x)*x-0.15*x;
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
double min(double* x1,double* x2)
{
	if(*x1 < *x2)
		return *x1;
	else
		return *x2;
}
double eval(double* a,double x)
{
	int i;
	double value=0.0;
	for(i=0;i<mo;i++)
		value += a[i]*fle(i,x);
	return value;
}
// the function of flux
double flux(double x)
{
	return -x;
}
void rk(double **u,double **hg,double *x,double *dx,int n)
{
	double v[n][mo];
	double temp1[n][mo];
	double temp2[n][mo];
	double temp3[n][mo];
	int i,j;

	for(i=0;i<mo;i++)
		for(j=0; j< n; j++){
			v[j][i]=u[j][i];
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
}
void res(double **u,double **hg,double *x,double *dx,int n)
{
	double h[n][mo];
	double w[n][mo];
	double un[n];
	double unn[n];
	double up[n];
	double bl[mo];
	double br[mo];
	int i,j,k;
	for(i=0;i<mo;i++)
	{
		bl[i]=fle(i,-0.5);
		br[i]=fle(i,0.5);
	}
	for(i=0;i<n;i++)
	{
		un[i]=eval((*(u+i)),0.5);
		up[i]=eval((*(u+i)),-0.5);
	}
	unn[0]=un[n-1];
	for(i=1;i<n;i++)
		unn[i]=un[i-1];
	for(i=0;i<n;i++)
		for(j=0;j<mo;j++)
			h[i][j]=un[i]*br[j]-unn[i]*bl[j];
	if(mp > 0)
			for(i=0;i<n;i++)
			{
				if(mp==1)
					h[i][1]=h[i][1]+(flux(up[i])+
							flux(eval((*(u+i)),0.0))*4.0+flux(un[i]))/6.0;
		/* 		else if(mp==2)
		 * 			for(k=1;k<mo;k++)
		 * 				h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*0.083333333+
		 * 						(eval((*(u+i)),-0.223606798)*fled(k,-0.223606798)+
		 * 	eval((*(u+i)),0.223606798)*fled(k,0.223606798))*0.416666667);
		 * 		else if(mp==3)
		 * 			for(k=1;k<mo;k++)
		 * 				h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*0.05+
		 * 						(flux(eval((*(u+i)),-0.458257569))*fled(k,-0.458257569)+
		 * 							flux(eval((*(u+i)),0.458257569))*fled(k,0.458257569))*0.272222222+
		 * 						(flux(eval((*(u+i)),0))*fle(k,0))*0.355555556);
		 */
				else
					for(k=1;k<mo;k++)
						h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*(1.0/30.0)+(flux(eval((*(u+i)),0.142615758))*fled(k,0.142615758)+flux(eval((*(u+i)),-0.142615758))*fled(k,-0.142615758))*0.277429189+(flux(eval((*(u+i)),0.382527662))*fled(k,0.382527662)+flux(eval((*(u+i)),-0.382527662))*fled(k,-0.382527662))*0.189237478);
			}
	for(i=0;i<n;i++)
		for(k=0;k<mo;k++)
			w[i][k]=ma[k]*h[i][k]/dx[i];
	for(i=0;i<n;i++)
	{
		un[i]=eval((*(w+i)),0.5);
		up[i]=eval((*(w+i)),-0.5);
	}
	unn[n-1]=up[0];
	for(i=0;i<n-1;i++)
		unn[i]=up[i+1];
	for(i=0;i<n;i++)
		for(j=0;j<mo;j++)
			h[i][j]=unn[i]*br[j]-up[i]*bl[j];
	if(mp > 0)
		for(i=0;i<n;i++)
		{
			if(mp==1)
				h[i][1]=h[i][1]+(flux(up[i])+
						flux(eval((*(w+i)),0.0))*4.0+flux(un[i]))/6.0;
	/* 		else if(mp==2)
	 * 			for(k=1;k<mo;k++)
	 * 				h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*0.083333333+
	 * 						(eval((*(w+i)),-0.223606798)*fled(k,-0.223606798)+
	 * 	eval((*(w+i)),0.223606798)*fled(k,0.223606798))*0.416666667);
	 * 		else if(mp==3)
	 * 			for(k=1;k<mo;k++)
	 * 				h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*0.05+(flux(eval((*(w+i)),-0.458257569))*fled(k,-0.458257569)+flux(eval((*(w+i)),0.458257569))*fled(k,0.458257569))*0.272222222+(flux(eval((*(w+i)),0))*fle(k,0))*0.355555556);
	 */
			else
				for(k=1;k<mo;k++)
					h[i][k]=h[i][k]+((flux(up[i])*fled(k,-0.5)+flux(un[i])*fled(k,0.5))*(1.0/30.0)+(flux(eval((*(w+i)),0.142615758))*fled(k,0.142615758)+flux(eval((*(w+i)),-0.142615758))*fled(k,-0.142615758))*0.277429189+(flux(eval((*(w+i)),0.382527662))*fled(k,0.382527662)+flux(eval((*(w+i)),-0.382527662))*fled(k,-0.382527662))*0.189237478);
		}
	for(i=0;i<n;i++)
		for(k=0;k<mo;k++)
			hg[i][k]=ma[k]*h[i][k]/dx[i];
}

double max(double* x1,double* x2)
{
	if(*x1 > *x2)
		return *x1;
	else
		return *x2;
}


