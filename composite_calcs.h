// ****************************************************************************
// file:        composite_calcs.h
// author:      Zane Rothe
// brief:       Contains function declarations for composite_calcs.c
// Class:       CPE 2600-111
// Assignment:  Lab 13 (final project)
// Date:        12/4/2024
// Compile:     makefile
// ****************************************************************************
#include <pthread.h>

typedef struct // data to send to calculation
{
    double t; //single layer thickness (in)
	double F1t; //psi
    double F1c; //psi
    double F2t; //psi
    double F2c; //psi
    double F6; //psi
    double BestSF; //Current best Safety Factor
    double N_M[6]; //Line loads and line moments
    double Q[9]; // Lamina Stiffness matrix
    int layers; // number of layers
    int iters; //current iteration number
    int iters_total; // total iterations
    int iters_perc; //rounded iteration number
    pthread_mutex_t mymutex; //locks common variables
    double BestSeq[]; //Current best stacking sequence (keep size flexible)
} data;


int nested_for(double init, double max, double step, double angles[], int rem, int current,data* mydata);
void halpin_tsai(double Pf[],double Pm[],double Pc[],double Vf);
void Qmat(double Pc[],double Q[]);
void Tmat(double theta,double T[]);
void Tstarmat(double theta,double Tstar[]);
double Det33(double mat[]);
void Invmat33(double mat[],double matinv[]);
void Multmat33x33(double a[],double b[],double c[]);
void ABDmats(double Q[],double t,double stack_seq[],int layers,double A[],double B[],double D[]);
void abcdmats(double A[],double B[],double C[],double D[],double result[]);
void Multmat66x61(double a[],double b[],double c[]);
void Multmat33x31(double a[],double b[],double c[]);
double tsai_hill(double S_12[],double F1t,double F1c,double F2t,double F2c,double F6);
double get_F1t(double Vf,double Ff1t,double Ef1,double Fmt,double Em1);
double get_F1c(double Vf,double Ef1,double Em1,double vf12,double vm12,double Gm12,double Fmt);
double get_K(double Vf,double Pm,double Pf);