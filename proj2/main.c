//
//  main.c
//  impsolve
//
//  Created by siqilu12345 on 16/11/25.
//  Copyright (c) 2016年 siqilu12345. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#define nprocess 5
#define total 1561
#define N 21
#define deltax 0.1
#define deltat 0.1
#define lambda 10
#define m 312
#define M 313
struct sparsematrix //defint the data structure of sparse matrix
{
    int coln[7];
    double val[7];
};
int cond(int i,int j,int k) /* judge if in the area*/
{
    if ((abs(i-(N-1)/2)+abs(j-(N-1)/2)+abs(k-(N-1)/2))<=(N-1)/2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void getrealxyz(int (*realxyz)[3],int n,int myrank) //get i,j,k from the number
{
    int realn=myrank*total/nprocess+n;
    int cout=0;
    int i,j,k;
    for (i=0; i<=N-1; i++)
    {
        for (j=0; j<=N-1; j++)
        {
            for (k=0; k<=N-1; k++)
            {
                if (cond(i,j,k)==1)
                {
                    if (cout==realn)
                    {
                        realxyz[n][0]=i;
                        realxyz[n][1]=j;
                        realxyz[n][2]=k;
                        cout++;
                    }
                    else
                    {
                        cout++;
                    }
                }
            }
        }
    }
}
void getresrealxyz(int (*resrealxyz)[3],int n,int myrank)
{
    int realn=myrank*total/nprocess+n+total/nprocess;
    int cout=0;
    int i,j,k;
    for (i=0; i<=N-1; i++)
    {
        for (j=0; j<=N-1; j++)
        {
            for (k=0; k<=N-1; k++)
            {
                if (cond(i,j,k)==1)
                {
                    if (cout==realn)
                    {
                        resrealxyz[n][0]=i;
                        resrealxyz[n][1]=j;
                        resrealxyz[n][2]=k;
                        cout++;
                    }
                    else
                    {
                        cout++;
                    }
                }
            }
        }
    }
}
int s1(int i,j,k) /*judge if at up surface*/
{
    if (k>=(N-1)/2&&(abs(i-(N-1)/2)+abs(j-(N-1)/2)+abs(k-(N-1)/2))==(N-1)/2)
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
    if (k<(N-1)/2&&(abs(i-(N-1)/2)+abs(j-(N-1)/2)+abs(k-(N-1)/2))==(N-1)/2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
double rhs(int a,int b,int c,int d) /*rhs of the equation*/
{
    double x,y,z,t;
    x=(a-(N-1)/2)*deltax;
    y=(b-(N-1)/2)*deltax;
    z=(c-(N-1)/2)*deltax;
    t=d*deltat;
    return 100*(sin(x)*sin(y)*sin(z)+3*t*sin(x)*sin(y)*sin(z));
}
int getxplusone(int x, int y, int z) /*get the number of (x+1,y,z)*/
{
    int i,j,k;
    int cout=0;
    for (i=x; i<=x+1; i++)
    {
        for (j=0; j<=N-1; j++)
        {
            for (k=0; k<=N-1; k++)
            {
                if (cond(i, j, k))
                {
                    if (i*(N+1)*(N+1)+j*(N+1)+k>x*(N+1)*(N+1)+y*(N+1)+z&&i*(N+1)*(N+1)+j*(N+1)+k<=(x+1)*(N+1)*(N+1)+y*(N+1)+z)
                    {
                        cout++;
                    }
                }
            }
        }
    }
    return cout;
}
int getyplusone(int x, int y, int z)
{
    int i,j,k;
    int cout=0;
    for (j=y; j<=y+1; j++)
    {
        for (k=0; k<=N-1; k++)
        {
            if (cond(i, j, k))
            {
                if (j*(N+1)+k>y*(N+1)+z&&j*(N+1)+k<=(y+1)*(N+1)+z)
                {
                    cout++;
                }
            }
        }
    }
    return cout;
}
int getzplusone(int x,y,z)
{
    return 1;
}
int getxminusone(int x, int y, int z)
{
    int i,j,k;
    int cout=0;
    for (i=x-1; i<=x; i++)
    {
        for (j=0; j<=N-1; j++)
        {
            for (k=0; k<=N-1; k++)
            {
                if (cond(i, j, k))
                {
                    if (i*(N+1)*(N+1)+j*(N+1)+k<x*(N+1)*(N+1)+y*(N+1)+z&&i*(N+1)*(N+1)+j*(N+1)+k>=(x-1)*(N+1)*(N+1)+y*(N+1)+z)
                    {
                        cout++;
                    }
                }
            }
        }
    }
    return -1*cout;
}
int getyminusone(int x, int y, int z)
{
    int i,j,k;
    int cout=0;
    for (j=y-1; j<=y; j++)
    {
        for (k=0; k<=N-1; k++)
        {
            if (cond(i, j, k))
            {
                if (j*(N+1)+k<y*(N+1)+z&&j*(N+1)+k>=(y-1)*(N+1)+z)
                {
                    cout++;
                }
            }
        }
    }
    return -1*cout;
}
int getzminusone(int x,y,z)
{
    return -1;
}
void giveb(double b[m],double resb[M-m],double u[m],double resu[m],int realxyz[m][3],int resrealxyz[M-m][3],int n) //give the value of b
{
    int myrank,i,j,k;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    for (i=0; i<=m-1; i++)
    {
        int x,y,z;
        x=realxyz[i][0];
        y=realxyz[i][1];
        z=realxyz[i][2];
        double a,bb,c,t;
        a=(x-(N-1)/2)*deltax;
        bb=(y-(N-1)/2)*deltax;
        c=(z-(N-1)/2)*deltax;
        t=n*deltat;
        if (s1(x, y, z))
        {
            b[i]=100*t*sin(a)*sin(bb)*sin(c);
        }
        else if (s2(x, y, z))
        {
            if (x>(N-1)/2&&y>(N-1)/2)
            {
                b[i]=100*(sqrt(3)*t*cos(a)*sin(bb)*sin(c)+sqrt(3)*t*sin(a)*cos(bb)*sin(c)-sqrt(3)*t*sin(a)*sin(bb)*cos(c));
            }
            if (x>(N-1)/2&&y<(N-1)/2)
            {
                b[i]=100*(sqrt(3)*t*cos(a)*sin(bb)*sin(c)-sqrt(3)*t*sin(a)*cos(bb)*sin(c)-sqrt(3)*t*sin(a)*sin(bb)*cos(c));
            }
            if (x<(N-1)/2&&y>(N-1)/2)
            {
                b[i]=100*(-sqrt(3)*t*cos(a)*sin(bb)*sin(c)+sqrt(3)*t*sin(a)*cos(bb)*sin(c)-sqrt(3)*t*sin(a)*sin(bb)*cos(c));
            }
            if (x<(N-1)/2&&y<(N-1)/2)
            {
                b[i]=100*(-sqrt(3)*t*cos(a)*sin(bb)*sin(c)-sqrt(3)*t*sin(a)*cos(bb)*sin(c)-sqrt(3)*t*sin(a)*sin(bb)*cos(c));
            }
            if (x==(N-1)/2&&y>(N-1)/2)
            {
                b[i]=100*(sqrt(2)*t*sin(a)*cos(bb)*sin(c)-sqrt(2)*t*sin(a)*sin(bb)*cos(c));
            }
            if (x==(N-1)/2&&y<(N-1)/2)
            {
                b[i]=100*(-sqrt(2)*t*sin(a)*cos(bb)*sin(c)-sqrt(2)*t*sin(a)*sin(bb)*cos(c));
            }
            if (x>(N-1)/2&&y==(N-1)/2)
            {
                b[i]=100*(sqrt(2)*t*cos(a)*sin(bb)*sin(c)-sqrt(2)*t*sin(a)*sin(bb)*cos(c));
            }
            if (x<(N-1)/2&&y==(N-1)/2)
            {
                b[i]=100*(-sqrt(2)*t*cos(a)*sin(bb)*sin(c)-sqrt(2)*t*sin(a)*sin(bb)*cos(c));
            }
            if (x==(N-1)/2&&y==(N-1)/2)
            {
                b[i]=-100*t*sin(a)*sin(bb)*cos(c);
            }
        }
        else
        {
            b[i]=rhs(x, y, z, n)*deltat+u[i];
        }
    }
    if (myrank==nprocess-1)
    {
        for (i=0; i<=M-m-1; i++)
        {
            int x,y,z;
            x=resrealxyz[i][0];
            y=resrealxyz[i][1];
            z=resrealxyz[i][2];
            double a,bb,c,t;
            a=(x-(N-1)/2)*deltax;
            bb=(y-(N-1)/2)*deltax;
            c=(z-(N-1)/2)*deltax;
            t=n*deltat;
            if (s1(x, y, z))
            {
                resb[i]=100*t*sin(a)*sin(bb)*sin(c);
            }
            else if (s2(x, y, z))
            {
                if (x>(N-1)/2&&y>(N-1)/2)
                {
                    resb[i]=100*(sqrt(3)*t*cos(a)*sin(bb)*sin(c)+sqrt(3)*t*sin(a)*cos(bb)*sin(c)-sqrt(3)*t*sin(a)*sin(bb)*cos(c));
                }
                if (x>(N-1)/2&&y<(N-1)/2)
                {
                    resb[i]=100*(sqrt(3)*t*cos(a)*sin(bb)*sin(c)-sqrt(3)*t*sin(a)*cos(bb)*sin(c)-sqrt(3)*t*sin(a)*sin(bb)*cos(c));
                }
                if (x<(N-1)/2&&y>(N-1)/2)
                {
                    resb[i]=100*(-sqrt(3)*t*cos(a)*sin(bb)*sin(c)+sqrt(3)*t*sin(a)*cos(bb)*sin(c)-sqrt(3)*t*sin(a)*sin(bb)*cos(c));
                }
                if (x<(N-1)/2&&y<(N-1)/2)
                {
                    resb[i]=100*(-sqrt(3)*t*cos(a)*sin(bb)*sin(c)-sqrt(3)*t*sin(a)*cos(bb)*sin(c)-sqrt(3)*t*sin(a)*sin(bb)*cos(c));
                }
                if (x==(N-1)/2&&y>(N-1)/2)
                {
                    resb[i]=100*(sqrt(2)*t*sin(a)*cos(bb)*sin(c)-sqrt(2)*t*sin(a)*sin(bb)*cos(c));
                }
                if (x==(N-1)/2&&y<(N-1)/2)
                {
                    resb[i]=100*(-sqrt(2)*t*sin(a)*cos(bb)*sin(c)-sqrt(2)*t*sin(a)*sin(bb)*cos(c));
                }
                if (x>(N-1)/2&&y==(N-1)/2)
                {
                    resb[i]=100*(sqrt(2)*t*cos(a)*sin(bb)*sin(c)-sqrt(2)*t*sin(a)*sin(bb)*cos(c));
                }
                if (x<(N-1)/2&&y==(N-1)/2)
                {
                    resb[i]=100*(-sqrt(2)*t*cos(a)*sin(bb)*sin(c)-sqrt(2)*t*sin(a)*sin(bb)*cos(c));
                }
                if (x==(N-1)/2&&y==(N-1)/2)
                {
                    resb[i]=-100*t*sin(a)*sin(bb)*cos(c);
                }
            }
            else
            {
                resb[i]=rhs(x, y, z, n)*deltat+resu[i];
            }
        }

    }
}
void jacobisolve(struct sparsematrix A[m],struct sparsematrix resA[M-m],double b[m],double resb[M-m],double u[m],double resu[M-m])//use jacobi to solve the equation
{
    int myrank,numprocess,flag=0;
    int i,j,k;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocess);
    double uu[total];
    for (k=0; k<=5000; k++)
    {
        MPI_Allgather(u, m, MPI_DOUBLE, uu, m, MPI_DOUBLE, MPI_COMM_WORLD);
        for (i=total-(M-m); i<=total-1; i++)
        {
            uu[i]=resu[i-(total-M+m)];
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for (i=0; i<=m-1; i++)
        {
            u[i]=0.0;
            for (j=1; j<=6; j++)
            {
                if (A[i].coln[j]!=-1)
                {
                    u[i]=u[i]-A[i].val[j]*uu[A[i].coln[j]]/A[i].val[0];
                }
            }
            u[i]=u[i]+b[i]/A[i].val[0];
        }
        if (myrank==nprocess-1)
        {
            for (i=0; i<=M-m-1; i++)
            {
                resu[i]=0.0;
                for (j=1; j<=6; j++)
                {
                    if (resA[i].coln[j]!=-1)
                    {
                        resu[i]=resu[i]-resA[i].val[j]*uu[resA[i].coln[j]]/resA[i].val[0];
                    }
                }
                resu[i]=resu[i]+resb[i]/resA[i].val[0];
            }
        }
        MPI_Bcast(resu, M-m, MPI_DOUBLE, nprocess-1, MPI_COMM_WORLD);
        flag=1;
        for (i=0; i<=M-m-1; i++)
        {
            if (fabs(resu[i]-uu[total-(M-m)+i])>=0.001)
            {
                flag=0;
                MPI_Bcast(&flag, 1, MPI_INTEGER, myrank, MPI_COMM_WORLD);
                break;
            }
        }
        for (i=0; i<=m-1; i++)
        {
            if (fabs(u[i]-uu[myrank*m+i])>=0.001)
            {
                flag=0;
                MPI_Bcast(&flag, 1, MPI_INTEGER, myrank, MPI_COMM_WORLD);
                break;
            }
        }
    }
}
int main(int argc, char *argv[])
{
    int i,j,k,myrank,numprocess;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    struct sparsematrix A[m],resA[M-m];
    double u[m],resu[M-m];
    double f[m],resf[M-m];
    double b[m],resb[M-m];
    int realxyz[m][3],resrealxyz[M-m][3];
    for (i=0; i<=m-1; i++)
    {
        getrealxyz(realxyz, i, myrank);
    }
    if (myrank==nprocess-1)
    {
        for (i=0; i<=M-m-1; i++)
        {
            getresrealxyz(resrealxyz,i,myrank);
        }
    }//give the value of realxyz
    for (i=0; i<=m-1; i++)
    {
        for (j=1; j<=6; j++)
        {
            A[i].coln[j]=-1;
        }
    }
    if (myrank==nprocess-1)
    {
        for (i=0; i<=M-m-1; i++)
        {
            for (j=1; j<=6; j++)
            {
                resA[i].coln[j]=-1;
            }
        }
    }
    for (i=0; i<=m-1; i++)
    {
        u[i]=0.0;
        f[i]=rhs(realxyz[i][0], realxyz[i][1], realxyz[i][2], 0);
        int t=myrank*m+i;
        if (s1(realxyz[i][0], realxyz[i][1], realxyz[i][2]))
        {
            A[i].coln[0]=t;
            A[i].val[0]=1;
        }
        else if (s2(realxyz[i][0], realxyz[i][1], realxyz[i][2]))
        {
            int x,y,z;
            x=realxyz[i][0];
            y=realxyz[i][1];
            z=realxyz[i][2];
            if (x>(N-1)/2&&y>(N-1)/2)
            {
                A[i].coln[0]=t;
                A[i].val[0]=3;
                A[i].coln[1]=t+getxminusone(x, y, z);
                A[i].val[1]=-1;
                A[i].coln[2]=t+getyminusone(x, y, z);
                A[i].val[2]=-1;
                A[i].coln[3]=t+getzplusone(x, y, z);
                A[i].val[3]=-1;
            }
            if (x>(N-1)/2&&y<(N-1)/2)
            {
                A[i].coln[0]=t;
                A[i].val[0]=3;
                A[i].coln[1]=t+getxminusone(x, y, z);
                A[i].val[1]=-1;
                A[i].coln[2]=t+getyplusone(x, y, z);
                A[i].val[2]=-1;
                A[i].coln[3]=t+getzplusone(x, y, z);
                A[i].val[3]=-1;
            }
            if (x<(N-1)/2&&y>(N-1)/2)
            {
                A[i].coln[0]=t;
                A[i].val[0]=3;
                A[i].coln[1]=t+getxplusone(x, y, z);
                A[i].val[1]=-1;
                A[i].coln[2]=t+getyminusone(x, y, z);
                A[i].val[2]=-1;
                A[i].coln[3]=t+getzplusone(x, y, z);
                A[i].val[3]=-1;
            }
            if (x<(N-1)/2&&y<(N-1)/2)
            {
                A[i].coln[0]=t;
                A[i].val[0]=3;
                A[i].coln[1]=t+getxplusone(x, y, z);
                A[i].val[1]=-1;
                A[i].coln[2]=t+getyplusone(x, y, z);
                A[i].val[2]=-1;
                A[i].coln[3]=t+getzplusone(x, y, z);
                A[i].val[3]=-1;
            }
            if (x==(N-1)/2&&y>(N-1)/2)
            {
                A[i].coln[0]=t;
                A[i].val[0]=2;
                A[i].coln[1]=t+getyminusone(x, y, z);
                A[i].val[1]=-1;
                A[i].coln[2]=t+getzplusone(x, y, z);
                A[i].val[2]=-1;
            }
            if (x==(N-1)/2&&y<(N-1)/2)
            {
                A[i].coln[0]=t;
                A[i].val[0]=2;
                A[i].coln[1]=t+getyplusone(x, y, z);
                A[i].val[1]=-1;
                A[i].coln[2]=t+getzplusone(x, y, z);
                A[i].val[2]=-1;
            }
            if (x>(N-1)/2&&y==(N-1)/2)
            {
                A[i].coln[0]=t;
                A[i].val[0]=2;
                A[i].coln[1]=t+getxminusone(x, y, z);
                A[i].val[1]=-1;
                A[i].coln[2]=t+getzplusone(x, y, z);
                A[i].val[2]=-1;
            }
            if (x<(N-1)/2&&y==(N-1)/2)
            {
                A[i].coln[0]=t;
                A[i].val[0]=2;
                A[i].coln[1]=t+getxplusone(x, y, z);
                A[i].val[1]=-1;
                A[i].coln[2]=t+getzplusone(x, y, z);
                A[i].val[2]=-1;
            }
            if (x==(N-1)/2&&y==(N-1)/2)
            {
                A[i].coln[0]=t;
                A[i].val[0]=1;
                A[i].coln[1]=t+getzplusone(x, y, z);
                A[i].val[1]=-1;
            }
        }
        else
        {
            int x,y,z;
            x=realxyz[i][0];
            y=realxyz[i][1];
            z=realxyz[i][2];
            A[i].coln[0]=t;
            A[i].val[0]=1+6*lambda;
            A[i].coln[1]=t+getxplusone(x, y, z);
            A[i].val[1]=-1*lambda;
            A[i].coln[2]=t+getxminusone(x, y, z);
            A[i].val[2]=-1*lambda;
            A[i].coln[3]=t+getyplusone(x, y, z);
            A[i].val[3]=-1*lambda;
            A[i].coln[4]=t+getyminusone(x, y, z);
            A[i].val[4]=-1*lambda;
            A[i].coln[5]=t+getzplusone(x, y, z);
            A[i].val[5]=-1*lambda;
            A[i].coln[6]=t+getzminusone(x, y, z);
            A[i].val[6]=-1*lambda;
        }
    }
    for (i=0; i<=M-m-1; i++)
    {
        resu[i]=0.0;
    }
    if (myrank==nprocess-1)
    {
        for (i=0; i<=M-m-1; i++)
        {
            int t=total-(M-m)+i;
            resf[i]=rhs(resrealxyz[i][0], resrealxyz[i][1], resrealxyz[i][2], 0);
            if (s1(resrealxyz[i][0], resrealxyz[i][1], resrealxyz[i][2]))
            {
                resA[i].coln[0]=i;
                resA[i].val[0]=1;
            }
            else if (s2(resrealxyz[i][0], resrealxyz[i][1], resrealxyz[i][2]))
            {
                int x,y,z;
                x=resrealxyz[i][0];
                y=resrealxyz[i][1];
                z=resrealxyz[i][2];
                if (x>(N-1)/2&&y>(N-1)/2)
                {
                    resA[i].coln[0]=t;
                    resA[i].val[0]=3;
                    resA[i].coln[1]=t+getxminusone(x, y, z);
                    resA[i].val[1]=-1;
                    resA[i].coln[2]=t+getyminusone(x, y, z);
                    resA[i].val[2]=-1;
                    resA[i].coln[3]=t+getzplusone(x, y, z);
                    resA[i].val[3]=-1;
                }
                if (x>(N-1)/2&&y<(N-1)/2)
                {
                    resA[i].coln[0]=t;
                    resA[i].val[0]=3;
                    resA[i].coln[1]=t+getxminusone(x, y, z);
                    resA[i].val[1]=-1;
                    resA[i].coln[2]=t+getyplusone(x, y, z);
                    resA[i].val[2]=-1;
                    resA[i].coln[3]=t+getzplusone(x, y, z);
                    resA[i].val[3]=-1;
                }
                if (x<(N-1)/2&&y>(N-1)/2)
                {
                    resA[i].coln[0]=t;
                    resA[i].val[0]=3;
                    resA[i].coln[1]=t+getxplusone(x, y, z);
                    resA[i].val[1]=-1;
                    resA[i].coln[2]=t+getyminusone(x, y, z);
                    resA[i].val[2]=-1;
                    resA[i].coln[3]=t+getzplusone(x, y, z);
                    resA[i].val[3]=-1;
                }
                if (x<(N-1)/2&&y<(N-1)/2)
                {
                    resA[i].coln[0]=t;
                    resA[i].val[0]=3.0;
                    resA[i].coln[1]=t+getxplusone(x, y, z);
                    resA[i].val[1]=-1.0;
                    resA[i].coln[2]=t+getyplusone(x, y, z);
                    resA[i].val[2]=-1.0;
                    resA[i].coln[3]=t+getzplusone(x, y, z);
                    resA[i].val[3]=-1.0;
                }
                if (x==(N-1)/2&&y>(N-1)/2)
                {
                    resA[i].coln[0]=t;
                    resA[i].val[0]=2.0;
                    resA[i].coln[1]=t+getyminusone(x, y, z);
                    resA[i].val[1]=-1.0;
                    resA[i].coln[2]=t+getzplusone(x, y, z);
                    resA[i].val[2]=-1.0;
                }
                if (x==(N-1)/2&&y<(N-1)/2)
                {
                    resA[i].coln[0]=t;
                    resA[i].val[0]=2.0;
                    resA[i].coln[1]=t+getyplusone(x, y, z);
                    resA[i].val[1]=-1.0;
                    resA[i].coln[2]=t+getzplusone(x, y, z);
                    resA[i].val[2]=-1.0;
                }
                if (x>(N-1)/2&&y==(N-1)/2)
                {
                    resA[i].coln[0]=t;
                    resA[i].val[0]=2.0;
                    resA[i].coln[1]=t+getxminusone(x, y, z);
                    resA[i].val[1]=-1.0;
                    resA[i].coln[2]=t+getzplusone(x, y, z);
                    resA[i].val[2]=-1.0;
                }
                if (x<(N-1)/2&&y==(N-1)/2)
                {
                    resA[i].coln[0]=t;
                    resA[i].val[0]=2.0;
                    resA[i].coln[1]=t+getxplusone(x, y, z);
                    resA[i].val[1]=-1.0;
                    resA[i].coln[2]=t+getzplusone(x, y, z);
                    resA[i].val[2]=-1.0;
                }
                if (x==(N-1)/2&&y==(N-1)/2)
                {
                    resA[i].coln[0]=t;
                    resA[i].val[0]=1.0;
                    resA[i].coln[1]=t+getzplusone(x, y, z);
                    resA[i].val[1]=-1.0;
                }
            }
            else
            {
                int x,y,z;
                x=resrealxyz[i][0];
                y=resrealxyz[i][1];
                z=resrealxyz[i][2];
                resA[i].coln[0]=t;
                resA[i].val[0]=1.0+6.0*lambda;
                resA[i].coln[1]=t+getxplusone(x, y, z);
                resA[i].val[1]=-1.0*lambda;
                resA[i].coln[2]=t+getxminusone(x, y, z);
                resA[i].val[2]=-1.0*lambda;
                resA[i].coln[3]=t+getyplusone(x, y, z);
                resA[i].val[3]=-1.0*lambda;
                resA[i].coln[4]=t+getyminusone(x, y, z);
                resA[i].val[4]=-1.0*lambda;
                resA[i].coln[5]=t+getzplusone(x, y, z);
                resA[i].val[5]=-1.0*lambda;
                resA[i].coln[6]=t+getzminusone(x, y, z);
                resA[i].val[6]=-1.0*lambda;
            }
        }
    }//give the value of A
    for (i=0; i<=20; i++)
    {
        giveb(b, resb, u, resu, realxyz, resrealxyz, i);
        jacobisolve(A, resA, b, resb, u, resu);
    } //move forward
    MPI_Finalize();
    return 0;
}
