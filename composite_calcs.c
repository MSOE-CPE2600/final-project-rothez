// ****************************************************************************
// file:        composite_calcs.c
// author:      Zane Rothe
// brief:       Contains functions for calculating composite properties
// Class:       CPE 2600-111
// Assignment:  Lab 13 (final project)
// Date:        12/4/2024
// Compile:     makefile
// ****************************************************************************
#include "composite_calcs.h"
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>

int nested_for(double init, double max, double step, double angles[],\
                                    int rem,int current,data* mydata)
// Uses nested loops to create for loop of custom depth
{
    current++; //increment current loop depth
    if (rem) //remaining depth
    {
        rem--;
        for(int i=init;i<max;i+=step)
        {
            angles[current]=i; // assemble stacking sequence
            nested_for(init,max,step,angles,rem,current,mydata); 
        }
    }
    else //at max depth in loops
    {
        // Do stuff here for each laminate combination
        double z[mydata->layers+1]; //Interface heights
        for(int i=0;i<=mydata->layers;i++)
        {
            z[i]=mydata->t*i-mydata->t*mydata->layers*0.5;
        }
        double A[9]={0};
        double B[9]={0};
        double D[9]={0};
        // Make A,B,D matrices for current stacking sequence
        ABDmats(mydata->Q,mydata->t,angles,mydata->layers,A,B,D);
        double abcd[36]={0};
        // Make a,b,c,d matrices
        abcdmats(A,B,B,D,abcd);
        double e_K[6];
        // find midplane strains and curvatures from stresses
        Multmat66x61(abcd,mydata->N_M,e_K);
        double layer_SF[mydata->layers];
        double e_xy[3];

        for(int i=0;i<mydata->layers;i++) //for each layer
        {
            for(int j=0;j<3;j++) // lamina plane strain
            {
                e_xy[j]=e_K[j]+(z[i]+z[i+1])/2*e_K[j+3];
            }
            double Tstar[9];
            Tstarmat(angles[i],Tstar);
            double e_12[3]; // Off-axis plane strain
            Multmat33x31(Tstar,e_xy,e_12);
            double S_12[3]; // Off-axis plane stress
            Multmat33x31(mydata->Q,e_12,S_12);
            layer_SF[i] = tsai_hill(S_12,mydata->F1t,mydata->F1c,\
                                mydata->F2t,mydata->F2c,mydata->F6);
        }
        double SF=layer_SF[0];
        for(int i=0;i<mydata->layers;i++)
        {
            if(layer_SF[i]<SF) // find lowest safety factor in layers
            {
                SF=layer_SF[i];
            }
        }

        pthread_mutex_lock(&mydata->mymutex); // dealing with shared variables
        mydata->iters++; //count iterations
        if((int)(100*(float)mydata->iters/mydata->iters_total)-1>\
                    mydata->iters_perc) //If additional 1% has finished
        {
            printf("\r\033[K");
            fflush(stdout);
            printf("%d%%",mydata->iters_perc); //show progress
            mydata->iters_perc++; //Percent of iterations done
            fflush(stdout);
        }
        if(SF>mydata->BestSF) // if new best found
        {
            mydata->BestSF=SF;
            for(int i=0;i<mydata->layers;i++)
            {
                mydata->BestSeq[i]=angles[i]; //save best sequence
            }
        }
        pthread_mutex_unlock(&mydata->mymutex);
    }
    return 0;
}

void halpin_tsai(double Pf[],double Pm[],double Pc[],double Vf)
// Calculates micromechanical properties of the composite from constituents
{
    double z[4]={1e9,2,1e9,1}; // Xi variable
    for(int i=0;i<4;i++)
    {
        double n=(Pf[i]/Pm[i]-1)/(Pf[i]/Pm[i]+z[i]); // Eta variable
        Pc[i]=(Pm[i]*(1+z[i]*n*Vf))/(1-n*Vf); // Composite property
    }
    return;
}

void Tmat(double theta,double T[])
// Calculates Transformation matrix T for given axis rotation angle theta
{
    double m=cos(theta*M_PI/180); //Theta in deg.
    double n=sin(theta*M_PI/180);
    T[0]=m*m;
    T[1]=n*n;
    T[2]=2*m*n;
    T[3]=n*n;
    T[4]=m*m;
    T[5]=-2*m*n;
    T[6]=-m*n;
    T[7]=m*n;
    T[8]=m*m-n*n;
    return;
}

void Tstarmat(double theta,double Tstar[])
// Calculates Transformation matrix T* for given axis rotation angle theta
{
    double m=cos(theta*M_PI/180); //Theta in deg.
    double n=sin(theta*M_PI/180);
    Tstar[0]=m*m;
    Tstar[1]=n*n;
    Tstar[2]=m*n;
    Tstar[3]=n*n;
    Tstar[4]=m*m;
    Tstar[5]=-m*n;
    Tstar[6]=-2*m*n;
    Tstar[7]=2*m*n;
    Tstar[8]=m*m-n*n;
    return;
}

void Qmat(double Pc[],double Q[])
// Caculates stiffness matrix Q for a lamina with properties Pc
{
    double delta=1-pow(Pc[2],2)*(Pc[1]/Pc[0]);
    Q[0]=Pc[0]/delta;
    Q[1]=Pc[2]*Pc[1]/delta;
    Q[2]=0;
    Q[3]=Q[1];
    Q[4]=Pc[1]/delta;
    Q[5]=0;
    Q[6]=0;
    Q[7]=0;
    Q[8]=Pc[3];
    return;
}

void Qbarmat(double Q[],double Qbar[],double T[],double Tstar[])
// Caculates off-axis stiffness matrix Qbar for lamina with stiffness matrix Q
{
    double Tinv[9];
    Invmat33(T,Tinv);
    double temp[9];
    Multmat33x33(Tinv,Q,temp);
    Multmat33x33(temp,Tstar,Qbar);
    return;
}

void Multmat33x33(double a[],double b[],double c[])
// Multiplies two 3x3 matrices
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            c[i+3*j]=0;
            for(int k=0;k<3;k++)
            {
                c[i+3*j]+=a[k+3*j]*b[i+3*k];
            }
        }
    }
    return;
}

void Multmat66x61(double a[],double b[],double c[])
// Multiplies 6x6 matrix with 6x1 matrix
{
    for(int i=0;i<6;i++)
    {
        c[i]=0;
        for(int k=0;k<6;k++)
        {
            c[i]+=a[k+6*i]*b[k];
        }
    }
    return;
}

void Multmat33x31(double a[],double b[],double c[])
// Multiplies 3x3 matrix with 3x1 matrix
{
    for(int i=0;i<3;i++)
    {
        c[i]=0;
        for(int k=0;k<3;k++)
        {
            c[i]+=a[k+3*i]*b[k];
        }
    }
    return;
}

double Det33(double mat[])
// Finds determinant of 3x3 matrix
{
    double det=(mat[0]*mat[4]*mat[8]+mat[1]*mat[5]*mat[6]+mat[2]*mat[3]*mat[7])\
    -(mat[6]*mat[4]*mat[2]+mat[7]*mat[5]*mat[0]+mat[8]*mat[3]*mat[1]);
    return det;
}

void Invmat33(double mat[],double matinv[])
// Inverts 3x3 matrix
{
    double det=Det33(mat);
    matinv[0]=(mat[4]*mat[8]-mat[5]*mat[7])/det;
    matinv[1]=-(mat[1]*mat[8]-mat[2]*mat[7])/det;
    matinv[2]=(mat[1]*mat[5]-mat[2]*mat[4])/det;
    matinv[3]=-(mat[3]*mat[8]-mat[5]*mat[6])/det;
    matinv[4]=(mat[0]*mat[8]-mat[2]*mat[6])/det;
    matinv[5]=-(mat[0]*mat[5]-mat[2]*mat[3])/det;
    matinv[6]=(mat[3]*mat[7]-mat[4]*mat[6])/det;
    matinv[7]=-(mat[0]*mat[7]-mat[1]*mat[6])/det;
    matinv[8]=(mat[0]*mat[4]-mat[1]*mat[3])/det;
    return;
}

void ABDmats(double Q[],double t,double stack_seq[],int layers,\
                                double A[],double B[],double D[])
// Calculates A,B,D matrices for laminate 
{
    double z[layers+1]; // height of each interface
    for(int i=0;i<=layers;i++)
    {
        z[i]=t*i-t*layers*0.5;
    }
    for(int i=0;i<layers;i++) // populate A,B,D matrices from each layer
    {
        double T[9];
        double Tstar[9];
        Tmat(stack_seq[i],T);
        Tstarmat(stack_seq[i],Tstar);
        double Qbar[9];
        Qbarmat(Q,Qbar,T,Tstar);
        for(int j=0;j<9;j++)
        {
            A[j]+=Qbar[j]*(z[i+1]-z[i]);
            B[j]+=(Qbar[j]/2)*(pow(z[i+1],2)-pow(z[i],2));
            D[j]+=(Qbar[j]/3)*(pow(z[i+1],3)-pow(z[i],3));
        }
    }
    return;
}

void abcdmats(double A[], double B[], double C[], double D[], double result[]){
    // create a,b,c,d matrices (6x6 inversion)
    // THE MATH IN THIS FUNCTION DEVELOPED WITH HELP FROM MICROSOFT COPILOT
    double invA[9], invD[9];
    Invmat33(A, invA);
    Invmat33(D, invD);

    // Compute intermediate matrices
    double invA_B[9], C_invA[9], C_invA_B[9], temp[9];
    Multmat33x33(invA, B, invA_B);
    Multmat33x33(C, invA, C_invA);
    Multmat33x33(C_invA, B, C_invA_B);

    // Compute the Schur complement
    for (int i = 0; i < 9; i++) {
        temp[i] = D[i] - C_invA_B[i];
    }

    double invSchur[9];
    Invmat33(temp, invSchur);

    // Compute top-left block: (A - B * invD * C)^-1
    double B_invD[9], B_invD_C[9], schurComplement[9], invSchurComplement[9];
    Multmat33x33(B, invD, B_invD);
    Multmat33x33(B_invD, C, B_invD_C);
    for (int i = 0; i < 9; i++) {
        schurComplement[i] = A[i] - B_invD_C[i];
    }
    Invmat33(schurComplement, invSchurComplement);

    // Fill the result matrix
    for (int i = 0; i < 36; i++) {
        result[i] = 0;
    }

    // Top-left 3x3 block
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i * 6 + j] = invSchurComplement[i * 3 + j];
        }
    }

    // Top-right 3x3 block (-invA * B * invSchur)
    double invA_B_invSchur[9];
    Multmat33x33(invA_B, invSchur, invA_B_invSchur);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i * 6 + (j + 3)] = -invA_B_invSchur[i * 3 + j];
        }
    }

    // Bottom-left 3x3 block (-invSchur * C * invA)
    double invSchur_C_invA[9];
    Multmat33x33(invSchur, C_invA, invSchur_C_invA);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[(i + 3) * 6 + j] = -invSchur_C_invA[i * 3 + j];
        }
    }

    // Bottom-right 3x3 block (invSchur)
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[(i + 3) * 6 + (j + 3)] = invSchur[i * 3 + j];
        }
    }

    return;
}

double tsai_hill(double S_12[],double F1t,double F1c,double F2t,double F2c,\
                                                                double F6)
// Uses Tsai-Hill failure criteria to determine safety factor
{
    double F1;
    double F2;
    if(S_12[0]>0)
    {
        F1=F1t;
    }
    else
    {
        F1=F1c;
    }

    if(S_12[1]>0)
    {
        F2=F2t;
    }
    else
    {
        F2=F2c;
    }

    double S_eff=pow(S_12[0]/F1,2)+pow(S_12[1]/F2,2)+\
        pow(S_12[2]/F6,2)-(S_12[0]*S_12[1]/pow(F1,2)); // Effective stress
    double SF=sqrt(1/S_eff);
    return SF;
}

double get_F1t(double Vf,double Ff1t,double Ef1,double Fmt,double Em1)
// Calculates Failure stress in longitudinal tension
{
    double F1t;
    double emt=Fmt/Em1; // matrix strain to failure
    double ef1t=Ff1t/Ef1; // fiber strain to failure
    if(emt>ef1t) // fibers more brittle
    {
        double Vf_crit=(Fmt-Em1*ef1t)/(Ff1t+Fmt-Em1*ef1t); // critical vol frac
        if(Vf>Vf_crit)
        {
            F1t=Ff1t*Vf+Em1*ef1t*(1-Vf);
        }
        else
        {
            F1t=Fmt*(1-Vf);
        }
    }
    else // matrix more brittle 
    {
        double Vf_crit=(Fmt)/(Ff1t+Fmt-Ef1*emt); // critical vol frac
        if(Vf>Vf_crit)
        {
            F1t=Ff1t*Vf;
        }
        else
        {
            F1t=Ef1*emt*Vf+Fmt*(1-Vf);
        }
    }
    return F1t;
}

double get_F1c(double Vf,double Ef1,double Em1,double vf12,\
                        double vm12,double Gm12,double Fmt)
// Calculates Failure stress in longitudinal compression
{
    double F1c;
    double OPB=2*Vf*sqrt((Em1*Ef1*Vf)/(3*(1-Vf))); // Out of phase buckling
    double IPB=Gm12/(1-Vf); // In phase buckling
    double TT=((Ef1*Vf+Em1*(1-Vf))*(1-cbrt(Vf))*(Fmt/Em1))/(vf12*Vf+vm12*(1-Vf));
            // Transverse tension
    if(OPB<IPB && OPB<TT)
    {
        F1c=OPB; // Out of phase buckling
    }
    else if (IPB<OPB && IPB<TT)
    {
        F1c=IPB; // In phase buckling
    }
    else
    {
        F1c=TT; // Transverse tension
    }
    return F1c;
}

double get_K(double Vf,double Pm,double Pf)
// Calcualtes stress concentraiton factor K for transverse tension/compression
{
    double K;
    K=(1-Vf*(1-Pm/Pf))/(1-sqrt(4*Vf/M_PI)*1-Pm/Pf);
    return K;
}