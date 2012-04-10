/*
 * =====================================================================================
 *
 *       Filename:  matrix.c
 *
 *    Description:  compute C = A x B, A is a matrix, B is a vector.The main process broadcast vector B to every slave process, then send matrix A by row to slave process, slave process compute the single row product B,and send the result to main process. the main process send A's row by turn,till end. once send final,the main process receive one result, the main process send a tag to slave process,the slave process will termination.the main process collect all data will be over.
 *
 *        Version:  1.0
 *        Created:  03/25/2012 08:38:46 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define  min(A,B)  ((A)<(B) ? (A) : (B))
#define  MAX_ROWS  1000
#define  MAX_COLS  1000
int main(int argc, char *argv[])
{
	int rows,cols;

	double a[MAX_ROWS][MAX_COLS],b[MAX_ROWS],c[MAX_ROWS];
	double buffer[MAX_ROWS],ans;

	int myid,master,numprocs;
	int i,j,numsent,numrcvd,sender;
	int anstype,row;
	double startwtime = 0.0, endwtime;
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	master = 0;
	rows = 100;
	cols = 100;
	int nsend = min(numprocs-1,rows);
	if(myid == master) {
    		startwtime = MPI_Wtime();
		//the main process initialize matrix A and vector B
		for(i = 0; i < cols; i++){
			b[i]=1.0;
			for(j = 0; j < rows; j++){
				a[j][i] = (double) i;
			}
		}
		numsent = 0;
		numrcvd = 0;
		//send b to the other slave process, by brodcast
		MPI_Bcast(b,rows,MPI_DOUBLE,master,MPI_COMM_WORLD);
		for(i = 0; i < nsend; i++) {
			for(j = 0; j < cols; j++)
				buffer[j] = a[i][j];
			// the (i+1)th process tackle with ith row,
			MPI_Send(buffer,cols,MPI_DOUBLE,i+1,i,MPI_COMM_WORLD);
			numsent++;
		}
		//all rows receive the result of the slave process computed
		for(i = 0; i < rows; i++){
			MPI_Recv(&ans,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			sender = status.MPI_SOURCE;
			anstype = status.MPI_TAG;
			c[anstype] = ans;
			// if it have more rows, compute continue
			if(numsent < rows){
				for(j=0; j < cols; j++)
					buffer[j] = a[numsent][j];
				MPI_Send(buffer,cols,MPI_DOUBLE,sender,numsent,MPI_COMM_WORLD);
				numsent++;
			}
			else // all rows has been tackled.
				//send a tag 0 message to slave process, let the slave process terminate
				MPI_Send((void *)NULL,0,MPI_DOUBLE,sender,9999,MPI_COMM_WORLD);
		}
		FILE *fp = fopen("result","a+");
		fprintf(fp,"The Result Vector C:\n");
		for(i = 0; i < rows; i++)
			fprintf(fp,"%f\n",c[i]);

	}
	else{
		// the below is the slave process
		// first, receive vector B
		MPI_Bcast(b,rows,MPI_DOUBLE,master,MPI_COMM_WORLD);
		// then, receive the one row in A,it is send by main process
		MPI_Recv(buffer,cols,MPI_DOUBLE,master,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		while(status.MPI_TAG != 9999) {
			// if the tag equal 0, exit
			row = status.MPI_TAG;
			ans = 0.0;
			for(i = 0; i < cols; i++)
				ans += buffer[i]*b[i];
			MPI_Send(&ans,1,MPI_DOUBLE,master,row,MPI_COMM_WORLD);
			MPI_Recv(buffer,cols,MPI_DOUBLE,master,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			// test
		}
	}
	if(myid == master)
        {
		endwtime = MPI_Wtime();
		printf("wall clock time = %f\n",endwtime-startwtime);
	}

	MPI_Finalize();
}


