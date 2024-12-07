// ****************************************************************************
// file:        composite.c
// author:      Zane Rothe
// brief:       Program to find optimal stacking sequence of composite material
// Args:        (1) number of layers, (2) file path to parameters
// Class:       CPE 2600-111
// Assignment:  Lab 13 (final project)
// Date:        12/4/2024
// Compile:     makefile
// Run example: $ ./composite 4 Comp_Input.csv 
// ****************************************************************************
#include "composite_calcs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <pthread.h>
#include <string.h>
#define MAX_LINE_LENGTH 50
#define MAX_LINES 24

int main(int argc, char* argv[])
{   
    // Read in parameters from file *******************************************
    printf("Welcome to Composite Compute!\n\nCalculating:\n");
    double input[MAX_LINES];
    int counter=0;
    FILE *file = fopen(argv[2], "r");
    if (!file) {
        perror("Unable to open file!");
        exit(1);
    }
    char line[MAX_LINE_LENGTH];
    while (fgets(line, MAX_LINE_LENGTH, file)) {
        strtok(line, ","); // Ignore label in first column
        char *number = strtok(NULL, ",");
        input[counter]=atof(number);
        counter++;
    }
    fclose(file);

    // Assign Parameter values to variables ***********************************
    // Laminate Dimensions (inches)
    double Lx=input[0];
    double Ly=input[1];
    double t=input[2];

    // Loads and moments (lbf), (in-lbf)
    double F_M[6]={input[3],input[4],input[5],input[6],input[7],input[8]};
    double N_M[6]={F_M[0]/Ly,F_M[1]/Lx,F_M[2]/Lx,F_M[3]/Ly,F_M[4]/Lx,F_M[5]/Lx};

    // Fiber properties (psi)
    double Ef1= input[9];
    double Ef2= input[10];
    double vf12= input[11];
    double Gf12= input[12];
    double Ff1t= input[13];
    double Pf[4]={Ef1,Ef2,vf12,Gf12};

    // Matrix properties (psi)
    double Em1= input[14];
    double Em2=Em1;
    double vm12= input[15];
    double Gm12=Em1/(2*(1+vm12));
    double Vf= input[16];
    double Fmt= input[17];
    double Fms= input[18];
    double Fmc= input[19];
    double Pm[4]={Em1,Em2,vm12,Gm12};

    // Calculate Lamina Failure Strengths (psi)
    double F1t=get_F1t(Vf,Ff1t,Ef1,Fmt,Em1);
    double F1c=get_F1c(Vf,Ef1,Em1,vf12,vm12,Gm12,Fmt);
    double F2t=Fmt/get_K(Vf,Em1,Ef2);
    double F2c=Fmc/get_K(Vf,Em1,Ef2);
    double F6 =Fms/get_K(Vf,Gm12,Gf12);

    //printf("%f %f %f %f %f\n",F1t,F1c,F2t,F2c,F6);

    int layers=atoi(argv[1]);

    // Micromechanical properties of composite
    double Pc[4];
    halpin_tsai(Pf,Pm,Pc,Vf);

    // Lamina Properties
    double Q[9];
    Qmat(Pc,Q);

    // Setup shared memory ****************************************************
    int shm_fd;
    pthread_mutex_t mymutex;
    pthread_mutexattr_t attr; // attribute
    pthread_mutexattr_init(&attr); // share mutex across processes
    pthread_mutexattr_setpshared(&attr, PTHREAD_PROCESS_SHARED);
    if (pthread_mutex_init(&mymutex,&attr))
    {
        perror("mutex_init");
        return 1;
    }
    //shared memory file descriptor
    shm_fd = shm_open("/myshm", O_CREAT | O_RDWR, 0666);
    if (shm_fd == -1)
    {
        perror("shm_open");
        return 1;
    }
    //set size of shared memory
    if (ftruncate(shm_fd, sizeof(data)+layers*sizeof(double)) == -1)
    {
        perror("ftruncate");
        return 1;
    }
    //data* mydata=malloc(sizeof(data)+layers*sizeof(double));
    data* mydata=mmap(0, sizeof(data)+layers*sizeof(double),\
             PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (mydata == MAP_FAILED)
    {
        perror("mmap");
        return 1;
    }

    // send data to struct for passing into loop ******************************
    mydata->F1t=F1t; //psi
    mydata->F1c=-F1c; //psi
    mydata->F2t=F2t; //psi
    mydata->F2c=-F2c; //psi
    mydata->F6=F6; //psi
    mydata->t=t; //in
    mydata->BestSF=0; //Current best Safety Factor
    mydata->layers=layers;
    mydata->mymutex=mymutex;
    
    for(int i=0;i<9;i++)
    {
        mydata->Q[i]=Q[i];
    }
    for(int i=0;i<6;i++)
    {
        mydata->N_M[i]=N_M[i];
    }

    // Begin multiple processes ***********************************************
    double ang_min = input[20];
    double ang_max = input[21];
    double resolution = input[22];
    int nprocs = 1; //find number of processes that satisfies resolution
    for (int i = 1; i <= ((ang_max-ang_min)/resolution); i++) 
    {
        if ((int)((ang_max-ang_min)/resolution) % i == 0) 
        {
            nprocs = i;
        }
        if (i>=(get_nprocs()-1))
        {
            break;
        }
    }
    double angles[layers];
    double ang_width=(ang_max-ang_min)/nprocs;
    double ang_min_div[nprocs];
    double ang_max_div[nprocs];

    mydata->iters = 0; // keeps track of progress
    mydata->iters_total = pow((ang_max-ang_min)/resolution,layers);
    
    for (int k=0;k<nprocs;k++)
    {
        ang_min_div[k]=ang_min+k*ang_width;
        ang_max_div[k]=ang_min_div[k]+ang_width;
        int pid=fork();
        if(pid==0)
        {
            for (double i=ang_min_div[k];i<ang_max_div[k];i+=resolution)
            {
                angles[0]=i;
                nested_for(ang_min,ang_max,resolution,angles,layers-1,0,mydata);
            } 
            exit(0); 
        }
    }
    for(int k=0;k<nprocs;k++)
    {
        wait(NULL);
    }
    // Results ****************************************************************
    printf("\r\033[K");
    fflush(stdout);
    printf("100%%\n");
    printf("\n\nBest Safety Factor: %.2f\n",mydata->BestSF);
    printf("Sequence: [");
    for(int i=0;i<layers;i++)
    {
        printf(" %.0f ",mydata->BestSeq[i]);
    }
    printf("]\n");
    // Cleanup ****************************************************************
    munmap(mydata, sizeof(data)+layers*sizeof(double));
    shm_unlink("/myshm");
    pthread_mutex_destroy(&mymutex);
    return 0;
}
/*
double ang_min = input[20];
    double ang_max = input[21];
    double resolution = input[22];


    mydata->iters = 0; // keeps track of progress
    mydata->iters_total = pow((ang_max-ang_min)/resolution+1,layers);
    
    double angles[layers];
    int nprocs = get_nprocs()-1; // use as many processors as are available
    double ang_width=(ang_max-ang_min)/nprocs;

    double ang_min_div[nprocs];
    double ang_max_div[nprocs];
    
    for (int k=0;k<nprocs;k++)
    {
        ang_min_div[k]=ang_min+k*ang_width;
        ang_max_div[k]=ang_min_div[k]+ang_width;
        int pid=fork();
        if(pid==0)
        {
            for (double i=ang_min_div[k];i<=ang_max_div[k];i+=resolution)
            {
                angles[0]=i;
                nested_for(ang_min,ang_max,resolution,angles,layers-1,0,mydata);
            } 
            exit(0); 
        }
    }
    for(int k=0;k<nprocs;k++)
    {
        wait(NULL);
    }
    */