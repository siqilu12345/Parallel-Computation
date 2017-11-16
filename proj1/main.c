//
//  main.c
//  para
//
//  Created by siqilu12345 on 16/10/31.
//  Copyright (c) 2016年 siqilu12345. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <malloc/malloc.h>
#define m 21 /*number of the points in one dimension*/
#define deltax 0.1 /*the step in one dimension, in space*/
#define deltat 0.1*0.1*0.05/*the step in time*/
#define lambda 0.05 /*the constant number of the mesh*/
double solution(double x,double y,double z,double t) /*the expression of the solution*/
{
    return 10000*t*sin(x)*sin(y)*sin(z);
}
int cond(int i,int j,int k) /* judge if in the area*/
{
    if ((abs(i-(m-1)/2)+abs(j-(m-1)/2)+abs(k-(m-1)/2))<=(m-1)/2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

double f(double x,double y,double z,double t) /*rhs of the equation*/
{
    return 10000*(sin(x)*sin(y)*sin(z)+3*t*sin(x)*sin(y)*sin(z));
}

double dirichlet(double x,double y,double z,double t) /*dirichlet condition*/
{
    return 10000*t*sin(x)*sin(y)*sin(z);
}

double neumannone(double x,double y,double z,double t) /* neumann condition on direction x*/
{
    return 10000*(sqrt(3)/3*t*cos(x)*sin(y)*sin(z)+sqrt(3)/3*t*sin(x)*cos(y)*sin(z)-sqrt(3)/3*t*sin(x)*sin(y)*cos(z));
}

double neumanntwo(double x,double y,double z,double t) /* neumann condition on direction x*/
{
    return 10000*(sqrt(3)/3*t*cos(x)*sin(y)*sin(z)-sqrt(3)/3*t*sin(x)*cos(y)*sin(z)-sqrt(3)/3*t*sin(x)*sin(y)*cos(z));
}

double neumannthree(double x,double y,double z,double t) /* neumann condition on direction x*/
{
    return 10000*(-1*sqrt(3)/3*t*cos(x)*sin(y)*sin(z)+sqrt(3)/3*t*sin(x)*cos(y)*sin(z)-sqrt(3)/3*t*sin(x)*sin(y)*cos(z));
}

double neumannfour(double x,double y,double z,double t) /* neumann condition on direction x*/
{
    return 10000*(-1*sqrt(3)/3*t*cos(x)*sin(y)*sin(z)-sqrt(3)/3*t*sin(x)*cos(y)*sin(z)-sqrt(3)/3*t*sin(x)*sin(y)*cos(z));
}
int s1(int i,j,k) /*judge if at up surface*/
{
    if (k>=(m-1)/2&&(abs(i-(m-1)/2)+abs(j-(m-1)/2)+abs(k-(m-1)/2))==(m-1)/2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int s2(int i,j,k) /*judge if at down surface*/
{
    if (k<(m-1)/2&&(abs(i-(m-1)/2)+abs(j-(m-1)/2)+abs(k-(m-1)/2))==(m-1)/2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void init(double *u[m][m],int h,int myrank,int nprocess)/*the initial status*/
{
    int i,j,k,z;
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=h-1; k++)
            {
                if (myrank>(nprocess-1)/2)
                {
                    z=(myrank-1)*h+k+1;
                }
                else if(myrank<(nprocess-1)/2)
                {
                    z=myrank*h+k;
                }
                else
                {
                    z=(m-1)/2;
                }/*get the real value z according to the rank of process*/
                if (cond(i, j, z))
                {
                    if (s1(i, j, z))/*on boundary*/
                    {
                        u[i][j][k]=0.0;
                    }
                    else if(s2(i, j, z))/*on boundary*/
                    {
                        u[i][j][k]=0.0;/*use the neumann boundary condition*/
                    }
                    else /*inside*/
                    {
                        u[i][j][k]=0.0;
                    }
                }
                else
                {
                    u[i][j][k]=0.0;
                }
            }
        }
    }
}
void nexttimeone(double *u[m][m], int h, double uprecv[m][m],double downrecv[m][m], int myrank,int nprocess,int n)
/* from the status right now to status after deltat time, without the down boundary*/
{
    int i,j,k,z;
    double *v[m][m];
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            v[i][j]=(double*)malloc((h+2)*sizeof(double));
        }
    }
    for (k=0; k<=h+1; k++)
    {
        for (i=0; i<=m-1; i++)
        {
            for (j=0; j<=m-1; j++)
            {
                if (k==0)
                {
                    v[i][j][k]=downrecv[i][j];
                }
                else if (k==h+1)
                {
                    v[i][j][k]=uprecv[i][j];
                }
                else
                {
                    v[i][j][k]=u[i][j][k-1];
                }
            }
        }/*record the status right now*/
    }
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=h-1; k++)
            {
                if (myrank>(nprocess-1)/2)
                {
                    z=(myrank-1)*h+k+1;
                }
                else if(myrank<(nprocess-1)/2)
                {
                    z=myrank*h+k;
                }
                else
                {
                    z=(m-1)/2;
                }/*get the real value z according to the rank of process*/
                if (cond(i, j, z))
                {
                    if (s1(i, j, z))/*on up boundary*/
                    {
                        u[i][j][k]=dirichlet((i-(m-1)/2)*deltax, (j-(m-1)/2)*deltax, (z-(m-1)/2)*deltax,n*deltat);
                    }
                    else if(s2(i, j, z)==0)/*inside*/
                    {
                        u[i][j][k]=(1-6*lambda)*v[i][j][k+1]+lambda*(v[i-1][j][k+1]+v[i+1][j][k+1]+v[i][j-1][k+1]+v[i][j+1][k+1]+v[i][j][k]+v[i][j][k+2])+deltat*f((i-(m-1)/2)*deltax,(j-(m-1)/2)*deltax,(z-(m-1)/2)*deltax,n*deltat);/*finite differential method*/
                    }
                }
            }
        }
    }
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            free(v[i][j]);
            v[i][j]=NULL;
        }
    }
}
void nexttimetwo(double *u[m][m], int h, double uprecv[m][m],double downrecv[m][m], int myrank,int nprocess,int n)
/* from the status right now to status after deltat time, use the neumann condition to calculate the down surface*/
{
    int i,j,k,z;
    double *v[m][m];
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            v[i][j]=(double*)malloc((h+2)*sizeof(double));
        }
    }
    for (k=0; k<=h+1; k++)
    {
        for (i=0; i<=m-1; i++)
        {
            for (j=0; j<=m-1; j++)
            {
                if (k==0)
                {
                    v[i][j][k]=downrecv[i][j];
                }
                else if (k==h+1)
                {
                    v[i][j][k]=uprecv[i][j];
                }
                else
                {
                    v[i][j][k]=u[i][j][k-1];
                }
            }
        }/*record the status right now*/
    }
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=h-1; k++)
            {
                if (myrank>(nprocess-1)/2)
                {
                    z=(myrank-1)*h+k+1;
                }
                else if(myrank<(nprocess-1)/2)
                {
                    z=myrank*h+k;
                }
                else
                {
                    z=(m-1)/2;
                }/*get the real value z according to the rank of process*/
                if (cond(i, j, z))
                {
                    if(s2(i, j, z))/*on down boundary*/
                    {
                        u[i][j][k]=0.0;
                        if (cond(i-1, j, z) && cond(i, j-1, z))
                        {
                            u[i][j][k]=(v[i-1][j][k+1]+v[i][j-1][k+1]+v[i][j][k+2])/3+deltax*sqrt(3)/3*neumannone((i-(m-1)/2)*deltax, (j-(m-1)/2)*deltax, (z-(m-1)/2)*deltax,n*deltat);
                        }
                        else if (cond(i-1, j, z) && cond(i, j+1, z))
                        {
                            u[i][j][k]=(v[i-1][j][k+1]+v[i][j+1][k+1]+v[i][j][k+2])/3+sqrt(3)/3*deltax*neumanntwo((i-(m-1)/2)*deltax, (j-(m-1)/2)*deltax, (z-(m-1)/2)*deltax,n*deltat);
                        }
                        else if (cond(i+1, j, z) && cond(i, j-1, z))
                        {
                            u[i][j][k]=(v[i+1][j][k+1]+v[i][j-1][k+1]+v[i][j][k+2])/3+deltax*sqrt(3)/3*neumannthree((i-(m-1)/2)*deltax, (j-(m-1)/2)*deltax, (z-(m-1)/2)*deltax,n*deltat);
                        }
                        else if (cond(i+1, j, z) && cond(i, j+1, z))
                        {
                            u[i][j][k]=(v[i+1][j][k+1]+v[i][j+1][k+1]+v[i][j][k+2])/3+deltax*sqrt(3)/3*neumannfour((i-(m-1)/2)*deltax, (j-(m-1)/2)*deltax, (z-(m-1)/2)*deltax,n*deltat);
                        }
                        /*use the neumann boundary condition*/
                    }
                }
            }
        }
    }
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            free(v[i][j]);
            v[i][j]=NULL;
        }
    }
}

int main(int argc, char *argv[])
{
    int nprocess, myrank;
    int i,j,k;
    int n=100; /*the number of repeats, indicates the time*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int h=(m-1)/(nprocess-1);/*processes do the similar things, except for the process rank (nprocess-1)/2*/
    double *u[m][m];
    if (myrank==(nprocess-1)/2)
    {
        h=1;
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=m-1;j++)
        {
            u[i][j]=(double*)malloc(h*sizeof(double));
        }
    }
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=h-1; k++)
            {
                u[i][j][k]=0.0;
            }
        }
    }
    double uprecv[m][m];
    double downrecv[m][m];
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            uprecv[i][j]=0.0;
            downrecv[i][j]=0.0;
        }
    }
    init(u, h, myrank, nprocess);
    for (k=1; k<=n; k++) /*solve the equation*/
    {
        if (myrank!=0 && myrank!=nprocess-1 && myrank!=(nprocess-1)/2)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    MPI_Sendrecv(&u[i][j][0], 1, MPI_DOUBLE, myrank-1, j, &downrecv[i][j], 1, MPI_DOUBLE, myrank-1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    MPI_Sendrecv(&u[i][j][h-1], 1, MPI_DOUBLE, myrank+1, j, &uprecv[i][j], 1, MPI_DOUBLE, myrank+1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                }
            }
        }
        else if (myrank==(nprocess-1)/2)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    MPI_Sendrecv(&u[i][j][0], 1, MPI_DOUBLE, myrank-1, j, &downrecv[i][j], 1, MPI_DOUBLE, myrank-1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    MPI_Sendrecv(&u[i][j][0], 1, MPI_DOUBLE, myrank+1, j, &uprecv[i][j], 1, MPI_DOUBLE, myrank+1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                }
            }
        }
        else if (myrank==0)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    MPI_Sendrecv(&u[i][j][h-1], 1, MPI_DOUBLE, myrank+1, j, &uprecv[i][j], 1, MPI_DOUBLE, myrank+1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    downrecv[i][j]=0.0;
                }
            }
        }
        else if (myrank==nprocess-1)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    MPI_Sendrecv(&u[i][j][0], 1, MPI_DOUBLE, myrank-1, j, &downrecv[i][j], 1, MPI_DOUBLE, myrank-1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    uprecv[i][j]=0.0;
                }
            }
        }
        nexttimeone(u, h, uprecv, downrecv, myrank, nprocess, k);
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank!=0 && myrank!=nprocess-1 && myrank!=(nprocess-1)/2)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    MPI_Sendrecv(&u[i][j][0], 1, MPI_DOUBLE, myrank-1, j, &downrecv[i][j], 1, MPI_DOUBLE, myrank-1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    MPI_Sendrecv(&u[i][j][h-1], 1, MPI_DOUBLE, myrank+1, j, &uprecv[i][j], 1, MPI_DOUBLE, myrank+1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                }
            }
        }
        else if (myrank==(nprocess-1)/2)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    MPI_Sendrecv(&u[i][j][0], 1, MPI_DOUBLE, myrank-1, j, &downrecv[i][j], 1, MPI_DOUBLE, myrank-1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    MPI_Sendrecv(&u[i][j][0], 1, MPI_DOUBLE, myrank+1, j, &uprecv[i][j], 1, MPI_DOUBLE, myrank+1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                }
            }
        }
        else if (myrank==0)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    MPI_Sendrecv(&u[i][j][h-1], 1, MPI_DOUBLE, myrank+1, j, &uprecv[i][j], 1, MPI_DOUBLE, myrank+1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    downrecv[i][j]=0.0;
                }
            }
        }
        else if (myrank==nprocess-1)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    MPI_Sendrecv(&u[i][j][0], 1, MPI_DOUBLE, myrank-1, j, &downrecv[i][j], 1, MPI_DOUBLE, myrank-1, j, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                    uprecv[i][j]=0.0;
                }
            }
        }
        nexttimetwo(u, h, uprecv, downrecv, myrank, nprocess, k);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    int z;
    for (i=0; i<=m-1; i++) /*print the result*/
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=h-1; k++)
            {
                if (myrank>(nprocess-1)/2)
                {
                    z=(myrank-1)*h+k+1;
                }
                else if(myrank<(nprocess-1)/2)
                {
                    z=myrank*h+k;
                }
                else
                {
                    z=(m-1)/2;
                }/*get the real value z according to the rank of process*/
                if (cond(i, j, z))
                {
                    printf("error u[%d][%d][%d]=%f from process %d\n",i,j,z,u[i][j][k]-solution((i-(m-1)/2)*deltax, (j-(m-1)/2)*deltax, (z-(m-1)/2)*deltax, n*deltat),myrank);
                }
            }
        }
    }
    MPI_Finalize();
    return 0;
}
