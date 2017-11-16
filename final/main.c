//
//  main.c
//  fftsolve
//
//  Created by siqilu12345 on 16/12/24.
//  Copyright (c) 2016年 siqilu12345. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define N 16
#define m 8
#define omega 0.923879533 - 0.382683432I
#define pi 3.14159265358979
void fftx(complex double u[m][m][m],int direction)/*fft on x direction*/
{
    int i,j,k;
    for (j=0; j<=m-1; j++)
    {
        for (k=0; k<=m-1; k++)
        {
            fftw_complex in[m],out[m];
            fftw_plan p;
            for (i=0; i<=m-1; i++)
            {
                in[i]=u[i][j][k];
            }
            if (direction==1)
            {
                p=fftw_plan_dft_1d(m, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            }
            else
            {
                p=fftw_plan_dft_1d(m, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
            }
            fftw_execute(p);
            for (i=0; i<=m-1; i++)
            {
                u[i][j][k]=out[i];
            }
            fftw_destroy_plan(p);
        }
    }
}
void exchangex(complex double even[m][m][m],complex double odd[m][m][m],int direction)/*use the data of even and odd to get the real fft on x direction*/
{
    complex double total[2*m][m][m];
    int i,j,k;
    complex double t;
    if (direction==1)
    {
        t=omega;
    }
    else
    {
        t=conj(omega);
    }
    for (j=0; j<=m-1; j++)
    {
        for (k=0; k<=m-1; k++)
        {
            for (i=0; i<=m-1; i++)
            {
                total[i][j][k]=even[i][j][k]+cpow(t, i)*odd[i][j][k];
                total[i+m][j][k]=even[i][j][k]-cpow(t, i)*odd[i][j][k];
            }
        }
    }
    for (j=0; j<=m-1; j++)
    {
        for (k=0; k<=m-1; k++)
        {
            for (i=0; i<=m-1; i++)
            {
                even[i][j][k]=total[2*i][j][k];
                odd[i][j][k]=total[2*i+1][j][k];
            }
        }
    }
}
void ffty(complex double u[m][m][m],int direction)/*fft on direction y*/
{
    int i,j,k;
    for (i=0; i<=m-1; i++)
    {
        for (k=0; k<=m-1; k++)
        {
            fftw_complex in[m],out[m];
            fftw_plan p;
            for (j=0; j<=m-1; j++)
            {
                in[j]=u[i][j][k];
            }
            if (direction==1)
            {
                p=fftw_plan_dft_1d(m, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            }
            else
            {
                p=fftw_plan_dft_1d(m, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
            }
            fftw_execute(p);
            for (j=0; j<=m-1; j++)
            {
                u[i][j][k]=out[j];
            }
            fftw_destroy_plan(p);
        }
    }
}
void exchangey(complex double even[m][m][m],complex double odd[m][m][m],int direction)/*use the data of even and odd to get the real fft on y direction*/
{
    complex double total[m][2*m][m];
    int i,j,k;
    complex double t;
    if (direction==1)
    {
        t=omega;
    }
    else
    {
        t=conj(omega);
    }
    for (i=0; i<=m-1; i++)
    {
        for (k=0; k<=m-1; k++)
        {
            for (j=0; j<=m-1; j++)
            {
                total[i][j][k]=even[i][j][k]+cpow(t, j)*odd[i][j][k];
                total[i][j+m][k]=even[i][j][k]-cpow(t, j)*odd[i][j][k];
            }
        }
    }
    for (i=0; i<=m-1; i++)
    {
        for (k=0; k<=m-1; k++)
        {
            for (j=0; j<=m-1; j++)
            {
                even[i][j][k]=total[i][2*j][k];
                odd[i][j][k]=total[i][2*j+1][k];
            }
        }
    }
}
void fftz(complex double u[m][m][m],int direction)/*fft on direction z*/
{
    int i,j,k;
    for (j=0; j<=m-1; j++)
    {
        for (i=0; i<=m-1; i++)
        {
            fftw_complex in[m],out[m];
            fftw_plan p;
            for (k=0; k<=m-1; k++)
            {
                in[k]=u[i][j][k];
            }
            if (direction==1)
            {
                p=fftw_plan_dft_1d(m, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            }
            else
            {
                p=fftw_plan_dft_1d(m, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
            }

            fftw_execute(p);
            for (k=0; k<=m-1; k++)
            {
                u[i][j][k]=out[k];
            }
            fftw_destroy_plan(p);
        }
    }
}
void exchangez(complex double even[m][m][m],complex double odd[m][m][m],int direction)/*use the data of even and odd to get the real fft on z direction*/
{
    complex double total[m][m][2*m];
    int i,j,k;
    complex double t;
    if (direction==1)
    {
        t=omega;
    }
    else
    {
        t=conj(omega);
    }
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=m-1; k++)
            {
                total[i][j][k]=even[i][j][k]+cpow(t, k)*odd[i][j][k];
                total[i][j][k+m]=even[i][j][k]-cpow(t, k)*odd[i][j][k];
            }
        }
    }
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=m-1; k++)
            {
                even[i][j][k]=total[i][j][2*k];
                odd[i][j][k]=total[i][j][2*k+1];
            }
        }
    }
}
int getrealz(int k,int myrank)/*the relationship of the real coordinate z with the process myrank and number k*/
{
    if (myrank%2==0)
    {
        return 2*k;
    }
    else
    {
        return 2*k+1;
    }
}
int getrealy(int j,int myrank)/*the relationship of the real coordinate y with the process myrank and number j*/
{
    if ((myrank/2)%2==0)
    {
        return 2*j;
    }
    else
    {
        return 2*j+1;
    }
}
int getrealx(int i,int myrank)/*the relationship of the real coordinate x with the process myrank and number i*/
{
    if (((myrank/2)/2)%2==0)
    {
        return 2*i;
    }
    else
    {
        return 2*i+1;
    }
}
void initf(complex double f[m][m][m],int myrank)/*function f*/
{
    int i,j,k;
    int x,y,z;
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=m-1; k++)
            {
                x=getrealx(i, myrank)+1;
                y=getrealy(j, myrank)+1;
                z=getrealz(k, myrank)+1;
                f[i][j][k]=3*sin(-pi+2*pi/N*x)*sin(-pi+2*pi/N*y)*sin(-pi+2*pi/N*z) +pow((sin(-pi+2*pi/N*x)*sin(-pi+2*pi/N*y)*sin(-pi+2*pi/N*z)),3);
            }
        }
    }
}
void initu(complex double u[m][m][m],int myrank)/*the initial of u*/
{
    int i,j,k;
    int x,y,z;
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=m-1; k++)
            {
                u[i][j][k]=0;
            }
        }
    }
}
double sumone(complex double u[m][m][m])/*sum of u*/
{
    int i,j,k;
    double a0=0.0;
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=m-1; k++)
            {
                a0=a0+creal(u[i][j][k]);
            }
        }
    }
    return a0;
}
double sumtwo(complex double u[m][m][m])/*sum of u^2*/
{
    int i,j,k;
    double a1=0.0;
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=m-1; k++)
            {
                a1=a1+creal(cpow(u[i][j][k],2));
            }
        }
    }
    return a1;
}
double sumthree(complex double u[m][m][m])/*sum of u^3*/
{
    int i,j,k;
    double a2=0.0;
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=m-1; k++)
            {
                a2=a2+creal(cpow(u[i][j][k],3));
            }
        }
    }
    return a2;
}
void initv(complex double v[m][m][m],complex double u[m][m][m])
{
    int i,j,k;
    for (i=0; i<=m-1; i++)
    {
        for (j=0; j<=m-1; j++)
        {
            for (k=0; k<=m-1; k++)
            {
                v[i][j][k]=cpow(u[i][j][k], 3);
            }
        }
    }
}
void fft3d(complex double f[m][m][m],complex double g[m][m][m],int direction)//3-dimension fft transform with mpi
{
    int nprocess,myrank,r,t;
    if (direction==1)
    {
        t=1;
    }
    else
    {
        t=0;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank<=7)//fft on direction x
    {
        fftx(f,t);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (r=0; r<=3; r++)//exchange data
    {
        if (myrank==r)
        {
            MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 0, MPI_COMM_WORLD);
        }
        if (myrank==r+4)
        {
            MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 1, MPI_COMM_WORLD);
        }
    }
    for (r=8; r<=11; r++)
    {
        if (myrank==r)
        {
            MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-8, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-4, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            exchangex(f, g, t);
        }
    }
    for (r=8; r<=11; r++)
    {
        if (myrank==r)
        {
            MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-8, 0, MPI_COMM_WORLD);
            MPI_Send(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-4, 1, MPI_COMM_WORLD);
        }
    }
    for (r=0; r<=3; r++)
    {
        if (myrank==r)
        {
            MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (myrank==r+4)
        {
            MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }//finish fft on direction x
    if (myrank<=7)//fft on direction y
    {
        ffty(f,t);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (r=0; r<=7; r++)//exchange data
    {
        if (r==0||r==1)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 0, MPI_COMM_WORLD);
            }
            if (myrank==r+2)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 1, MPI_COMM_WORLD);
            }
        }
        if (r==4||r==5)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+6, 0, MPI_COMM_WORLD);
            }
            if (myrank==r+2)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+6, 1, MPI_COMM_WORLD);
            }
        }
    }
    for (r=8; r<=11; r++)
    {
        if (r==8||r==9)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-8, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-6, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                exchangey(f, g, t);
            }
        }
        else
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-6, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-4, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                exchangey(f, g, t);
            }
        }
    }
    for (r=8; r<=11; r++)
    {
        if (r==8||r==9)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-8, 0, MPI_COMM_WORLD);
                MPI_Send(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-6, 1, MPI_COMM_WORLD);
            }
        }
        else
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-6, 0, MPI_COMM_WORLD);
                MPI_Send(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-4, 1, MPI_COMM_WORLD);
            }
            
        }
    }
    for (r=0; r<=7; r++)
    {
        if (r==0||r==1)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (myrank==r+2)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }
        if (r==4||r==5)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+6, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (myrank==r+2)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+6, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }
    }//finish fft on direction y
    if (myrank<=7)//fft on direction z
    {
        fftz(f,t);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (r=0; r<=7; r++)//exchange data
    {
        if (r==0)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 0, MPI_COMM_WORLD);
            }
            if (myrank==r+1)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 1, MPI_COMM_WORLD);
            }
        }
        if (r==2)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+7, 0, MPI_COMM_WORLD);
            }
            if (myrank==r+1)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+7, 1, MPI_COMM_WORLD);
            }
        }
        if (r==4)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+6, 0, MPI_COMM_WORLD);
            }
            if (myrank==r+1)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+6, 1, MPI_COMM_WORLD);
            }
        }
        if (r==6)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+5, 0, MPI_COMM_WORLD);
            }
            if (myrank==r+1)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+5, 1, MPI_COMM_WORLD);
            }
        }
    }
    for (r=8; r<=11; r++)
    {
        if (r==8)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-8, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-7, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                exchangez(f, g, t);
            }
        }
        if (r==9)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-7, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-6, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                exchangez(f, g, t);
            }
        }
        if (r==10)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-6, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-5, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                exchangez(f, g, t);
            }
        }
        if (r==11)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-5, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-4, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                exchangez(f, g, t);
            }
        }
    }
    for (r=8; r<=11; r++)
    {
        if (r==8)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-8, 0, MPI_COMM_WORLD);
                MPI_Send(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-7, 1, MPI_COMM_WORLD);
            }
        }
        if (r==9)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-7, 0, MPI_COMM_WORLD);
                MPI_Send(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-6, 1, MPI_COMM_WORLD);
            }
        }
        if (r==10)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-6, 0, MPI_COMM_WORLD);
                MPI_Send(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-5, 1, MPI_COMM_WORLD);
            }
        }
        if (r==11)
        {
            if (myrank==r)
            {
                MPI_Send(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r-5, 0, MPI_COMM_WORLD);
                MPI_Send(g, m*m*m, MPI_C_DOUBLE_COMPLEX, r-4, 1, MPI_COMM_WORLD);
            }
        }
    }
    for (r=0; r<=7; r++)
    {
        if (r==0)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (myrank==r+1)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+8, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }
        if (r==2)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+7, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (myrank==r+1)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+7, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }
        if (r==4)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+6, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (myrank==r+1)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+6, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }
        if (r==6)
        {
            if (myrank==r)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+5, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (myrank==r+1)
            {
                MPI_Recv(f, m*m*m, MPI_C_DOUBLE_COMPLEX, r+5, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }
    }//finish fft on direction z
}
double solution(double a3,double a2,double a1,double a0)//solve cubic equation
{
    complex double a,b,c,d;
    a=a3+0.0*I;
    b=a2+0.0*I;
    c=a1+0.0*I;
    d=a0+0.0*I;
    complex double delta,delta0,delta1,cc;
    complex double kasai;
    kasai=-0.5+sqrt(3)/2*I;
    complex double x[3];
    delta=18*a*b*c*d-4*cpow(b, 3)*d+cpow(b, 2)*cpow(c, 2)-4*a*cpow(c, 3)-27*cpow(a, 2)*cpow(d, 2);
    delta0=cpow(b, 2)-3*a*c;
    delta1=2*cpow(b, 3)-9*a*b*c+27*cpow(a, 2)*d;
    cc=cpow(delta1, 2)-4*cpow(delta0,3);
    cc=cpow((delta1+cpow(cc, 0.5))*0.5, 1.0/3);
    double s;
    int i;
    if (cabs(cc)<=0.000001)
    {
        s=creal(-b/3/a);
    }
    else
    {
        for (i=0; i<=2; i++)
        {
            x[i]=-1/(3*a)*(b+cpow(kasai, i)*cc+delta0/(cpow(kasai, i)*cc));
        }
        int flag=0;
        for (i=0; i<=2; i++)
        {
            if (fabs(cimag(x[i]))<=0.00001)
            {
                if (flag==0)
                {
                    s=creal(x[i]);
                    flag=1;
                }
                else
                {
                    if (cabs(x[i])<fabs(s))
                    {
                        s=creal(x[i]);
                    }
                }
            }
        }
    }
    return s;
}
int main(int argc, char *argv[])
{
    int nprocess,myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int i,j,k,t;
    int x,y,z,xx,yy,zz;
    double c0,a,b,c,a0,a1,a2,a3;
    complex double u[m][m][m],v[m][m][m],f[m][m][m],g[m][m][m];
    if (myrank<=7)//give f
    {
        initf(f, myrank);
        a=sumone(f);
    }
    MPI_Allreduce(&a, &a0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//get the sum of f
    if (myrank<=7)
    {
        initu(u, myrank);
    }
    fft3d(f, g, 1);//fft of f
    for (t=0; t<=30; t++)//solve the equation
    {
        if (myrank<=7)
        {
            initv(v, u);//v=u^3
        }
        fft3d(v, g, 1);//fft of u
        if (myrank<=7)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    for (k=0; k<=m-1; k++)
                    {
                        x=getrealx(i, myrank);
                        y=getrealy(j, myrank);
                        z=getrealz(k, myrank);
                        if (x==0&&y==0&&z==0)
                        {
                            u[i][j][k]=0.0;
                        }
                        else
                        {
                            if (x>=m+1)
                            {
                                xx=2*m-x;
                            }
                            else
                            {
                                xx=x;
                            }
                            if (y>=m+1)
                            {
                                yy=2*m-y;
                            }
                            else
                            {
                                yy=y;
                            }
                            if (z>=m+1)
                            {
                                zz=2*m-z;
                            }
                            else
                            {
                                zz=z;
                            }
                            u[i][j][k]=(f[i][j][k]-v[i][j][k])/(xx*xx+yy*yy+zz*zz)/N/(2*pi)/(2*pi);
                        }
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        fft3d(u, g, 0);//ifft of u
        if (myrank<=7)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    for (k=0; k<=m-1; k++)
                    {
                        u[i][j][k]=u[i][j][k]/N/N/N;
                    }
                }
            }
        }
        if (myrank<=7)
        {
            a=sumone(u);
            b=sumtwo(u);
            c=sumthree(u);
        }
        MPI_Allreduce(&a, &a3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&b, &a2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&c, &a1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank==8)
        {
            c0=solution(N*N*N+0.0,3*a3,3*a2,a1-a0);
        }//solve the cubic equation
        MPI_Bcast(&c0, 1, MPI_DOUBLE, 8, MPI_COMM_WORLD);
        if (myrank<=7)
        {
            for (i=0; i<=m-1; i++)
            {
                for (j=0; j<=m-1; j++)
                {
                    for (k=0; k<=m-1; k++)
                    {
                        u[i][j][k]=u[i][j][k]+c0;
                    }
                }
            }
        }
    }
    MPI_Finalize();
    return 0;
}
