/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <mpi.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
//#include <fortran.h>
#include <stdio.h>
#include <stdlib.h>
#include "seism_io_consts.h"
#include "seism_io.h"

int X,Y,Z;

//copy from pmcl3d.c
const double   micro = 1.0e-6;
double gethrtime()
{
    struct timeval TV;
    int RC = gettimeofday(&TV,NULL);

    if (RC == -1){
       printf("Bad call to gettimeofday\n");
       return(-1);
    }

    return ( ((double)TV.tv_sec ) + micro * ((double)  TV.tv_usec));
}

void cleanChunk(float *u)
{
    memset((void*)u, 0, X*Y*Z*sizeof(float));
}

void printChunk(float u[][Y][Z],int nxt, int nyt, int nzt,int ghostx, int ghosty, int ghostz, int rank)
{
    char info[1000000]="\0", tmp[200]="\0";

    sprintf(tmp,"%d ) Read Grid is:\n",rank);
    strcpy(info,tmp);
    for(int kk=ghostz; kk<nzt+ghostz*2-ghostz; kk++)
    for(int jj=ghosty; jj<nyt+ghosty*2-ghosty; jj++)
    {
       sprintf(tmp,"%3d )%3dx%3d=",rank, kk,jj);
       strcat(info,tmp);

       for(int ii=ghostx; ii<nxt+ghostx*2-ghostx; ii++)
       {
          sprintf(tmp," %5.1f ",u[ii][jj][kk]);
          strcat(info,tmp);
       }
       sprintf(tmp,"\n");
       if(nxt< 5) strcat(info,tmp);
    }
    printf("%s", info);
}

int test(MPI_Comm Comm, int rank , int size, int type, int cores, int chunksize)
{
    MPI_Comm MC1;
    int rankmc1;
    int err;

    long   start,end;
    double t1,t2,dtime;
    const int maxdim = 3;
    int coords[]={0,0,0};

    char seism_method[16];
    char sxgro[32];
    char sxgro2[32];

    //select type mode
    //int type = 3;
    if(type == 0){ /*MPIIO*/
        sprintf(seism_method, "mpiio");
        sprintf(sxgro, "output_sfc/CSX96PS_%d",size);
        sprintf(sxgro2, "output_vlm/CSX96PV_%d",size);
    }else if(type == 1){ /*PHDF5*/
        sprintf(seism_method, "phdf5");
        sprintf(sxgro, "output_sfc/CSX96PS_%d.h5",size);
        sprintf(sxgro2, "output_vlm/CSX96PV_%d.h5",size);
    }else if(type == 2){ /*PNETCDF*/
        sprintf(seism_method, "pnetcdf");
        sprintf(sxgro, "output_sfc/CSX96PS_%d.nc",size);
        sprintf(sxgro2, "output_vlm/CSX96PV_%d.nc",size);
    }else if(type == 3){ /*ADIOS*/
        sprintf(seism_method, "adios");
        sprintf(sxgro, "output_sfc/CSX96PS_%d.bp",size);
        sprintf(sxgro2, "output_vlm/CSX96PV_%d.bp",size);
    }else{
        printf("Invalid type.\n");
        return -2;
    }

    const int TIMES = 4;

    int nx,  ny,  nz;
    int nxt, nyt, nzt;
    int npx, npy, npz;

  if(cores== 1){
        ///1 cores
        npx=1; npy=1; npz=1;
    }else if(cores== 2){
        ///2 cores
        npx=1; npy=2; npz=1;
    }else if(cores== 4){
        ///4 cores
        npx=2; npy=2; npz=1;
    }else if(cores== 6){
        ///6 cores
        npx=3; npy=2; npz=1;
    }else if(cores== 8){
        ///8 cores
        npx=2; npy=2; npz=2;
    }else if(cores== 12){
        ///12 cores
        npx=3; npy=2; npz=2;
    }else if(cores== 16){
        ///16 cores
        npx=4; npy=2; npz=2;
    }else if(cores== 24){
        ///24 cores
       npx=3; npy=4; npz=2;
    }else if(cores== 36){
        ///36 cores
       npx=3; npy=4; npz=3;
    }else if(cores== 48){
        ///48 cores
        npx=3; npy=4; npz=4;
    }else if(cores== 96){
        ///96 cores
        npx=6; npy=4; npz=4;
    }else if(cores== 576){
        ///576 cores
        npx=8; npy=18; npz=4;
    }else if(cores== 1728){
        ///1728 cores
        npx=24; npy=18; npz=4;
    }else if(cores== 1024){
        ///1K cores
        npx=4; npy=16; npz=16;
    }else if(cores== 1152){
        ///1152 cores
        npx=16; npy=18; npz=4;
    }else if(cores== 2048){
        ///2K cores
        npx=8; npy=16; npz=16;
    }else if(cores== 3072){
        ///3K cores
        npx=12; npy=16; npz=16;
    }else if(cores== 4096){
        ///4K cores
        npx=16; npy=16; npz=16;
    }else if(cores== 10240){
        ///10K cores
        npx=40; npy=16; npz=16;
    }else if(cores== 24000){
        ///24000 cores 72Node Comet
        npx=60; npy=25; npz=16;
    }else if(cores== 51200){
        ///51200 cores 72Node Bluewaters
        npx=64; npy=50; npz=16;
    }else{
        ///default
        npy=16; npz=4;
        npx= cores/npy/npz;
    }

    if(chunksize==2){
        //small
        nxt=2; nyt=3; nzt=4;
    }else if(chunksize==99){
        //for flexiabe test
        nxt=120; nyt=130; nzt=3072;
    }else if(chunksize==8){
        //8MB
        nxt=128; nyt=128; nzt=128;
    }else if(chunksize==5){
        //5MB
        nxt=128; nyt=64; nzt=180;
    }else if(chunksize==112){
        //112MB
        nxt=120; nyt=80; nzt=3072;
    }else if(chunksize==200){
        //200MB
        nxt=160; nyt=160; nzt=2048;
    }else if(chunksize==400){
        //more biger
        nxt=320; nyt=160; nzt=2048;
    }else{
        //default
        nxt=2; nyt=3; nzt=4;
    }

    nx =npx * nxt;
    ny =npy * nyt;
    nz =npz * nzt;

    int ghostx=3, ghosty=3, ghostz=2;
    int dims[maxdim], periodic[maxdim];
    dims[0]=npx; dims[1]=npy; dims[2]=npz;
    periodic[0]=0; periodic[1]=0; periodic[2]=0;
    int reorder = 1;
    if(nx%nxt!= 0 || ny%nyt!=0 || nz%nzt!=0)
    {
       if(rank==0)
         printf("nx,ny,nz (%d,%d,%d) and nxt,nyt,nzt (%d,%d,%d) are not CORRECT.\nProgram will EXIT.\n", nx,ny,nz,nxt,nyt,nzt);
       //MPI_Barrier(Comm);
       //MPI_Finalize();
       return -1;
    }
    if(size != npx*npz*npy)
    {
       if(rank==0)
         printf("Size of the communicator (%d) is not equal size of the Cartesian topology (%d)\nProgram will EXIT.\n", size,npx*npy*npz);
       //MPI_Barrier(Comm);
       //MPI_Finalize();
       return -1;
    }
    else{
       if(rank==0)
         printf("Size:(%d), SEISM-IO TEST Cartesian Topology (%dx%dx%d). Method:%s\n", size,npx,npy,npz,seism_method);
    }

    MPI_Cart_create(Comm, maxdim, dims, periodic, reorder, &MC1);
    MPI_Cart_get(MC1, maxdim, dims, periodic, coords);
    MPI_Cart_rank(MC1, coords, &rankmc1);

    //if(rank==0) printf("%d) Cart coords=%d\t%d\t%d\n",rank,coords[0],coords[1],coords[2]);

    seism_init(&MC1,&rank,coords,(int *)&maxdim,
               &nx,&ny,&nz,&nxt,&nyt,&nzt,
               &ghostx,&ghosty,&ghostz,
               &npx,&npy,&npz,
               seism_method,&err);
    if(err < 0){
        printf("SEISM ERROR! Init failed!");
        MPI_Abort(Comm,-1);
        MPI_Finalize();
        return 1;
    }
    if(err > 0){
        printf("ERROR! %s is not defined in compilation, please check Makefile.\n", seism_method);
        MPI_Abort(Comm,-1);
        MPI_Finalize();
        return 1;
    }
    else{
        //if(rankmc1 ==0) printf("%d In main(), inited\n", rankmc1);
    }

    int nbgx=1, nedx=nx, nskpx=1,nbgy=1, nedy=ny, nskpy=1;
    int nbgz=1, nedz=1, nskpz=1;
    //nedz= nz-nedz+1;
    //nbgz= nz-nbgz+1;
    int seism_regGridID, seism_regGridID4_read;
    int seism_filex, seism_filex4_read;
    seism_createRegularGrid( &nbgx, &nedx, &nskpx, &nbgy, &nedy, &nskpy,  &nbgz,&nedz,
                             &nskpz, &seism_regGridID, &err);
    seism_createRegularGrid( &nbgx, &nedx, &nskpx, &nbgy, &nedy, &nskpy,  &nbgz,&nedz,
                             &nskpz, &seism_regGridID4_read, &err);


    int nbgx2=1, nedx2=nx, nskpx2=1,nbgy2=1, nedy2=ny, nskpy2=1;
    int nbgz2=1, nedz2=nz, nskpz2=1;
    //nedz2= nz-nedz2+1;
    //nbgz2= nz-nbgz2+1;
    int seism_regGridID2,seism_regGridID3_read;
    int seism_filex2, seism_filex3_read;
    seism_createRegularGrid( &nbgx2, &nedx2, &nskpx2, &nbgy2, &nedy2, &nskpy2,
                             &nbgz2,&nedz2, &nskpz2, &seism_regGridID2, &err);
    seism_createRegularGrid( &nbgx2, &nedx2, &nskpx2, &nbgy2, &nedy2, &nskpy2,
                             &nbgz2,&nedz2, &nskpz2, &seism_regGridID3_read, &err);

    if(rankmc1==0){
        //if(rankmc1 ==0) printf("%d In main(), created regulargrid 1,2,3.\n", rankmc1);
    }

    char fmode[]="w";
    char fdata[]="float";
    int write_step=1;
    seism_filex = -1111;

    time(&start);
    t1=MPI_Wtime();
    //if(rankmc1 ==0) printf("%d) 1 before open,%d, %lu\n",rank, seism_filex,strlen(sxgro));
    seism_file_open(sxgro, fmode, &write_step, fdata, &seism_regGridID, &seism_filex, &err);
    if(err !=0)
    {
        printf("SEISM ERROR! Open %s failed.", sxgro);
        MPI_Abort(MC1,-1);
        MPI_Finalize();
        return 1;
    }
    else{
        if(rankmc1 ==0) printf("%d) In main(), %s opened.\n", rankmc1, sxgro);
    }

    seism_filex2 = -1111;
    seism_file_open(sxgro2, fmode, &write_step, fdata, &seism_regGridID2, &seism_filex2, &err);

    if(err !=0)
    {
        printf("SEISM ERROR! Open %s failed.", sxgro2);
        MPI_Abort(MC1,-1);
        MPI_Finalize();
        return 1;
    }
    else{
        if(rankmc1 ==0) printf("%d) In main(), %s opened.\n", rankmc1, sxgro2);
    }

    t2=MPI_Wtime();
    MPI_Barrier(MC1);
    time(&end);
    t1 = t2-t1;
    MPI_Allreduce( &t1, &dtime, 1, MPI_DOUBLE, MPI_MAX, MC1);

    if(rankmc1 ==0) printf("%d) rank, Open Times: %ld(c native) %.1f(MPI_MAX) seconds.\n",rankmc1, end-start, dtime);

    //MPI_Barrier(Comm);
    //if(DEBUG>0 && rankmc1 ==0) printf("%d) In main(), before call seism_file_open.\n",rank);


    int sum;
    //initial chunck
    //float u1[nxt+ghostx*2][nyt+ghosty*2][nzt+ghostz*2];
    float *u1 = malloc((nxt+ghostx*2)*(nyt+ghosty*2)*(nzt+ghostz*2)*sizeof(float));
    if(u1 ==NULL)
    {
        if(rankmc1==0) printf("SEISM ERROR! Memory allocate failed.");
        MPI_Abort(Comm,-1);
        MPI_Finalize();
        return 1;
    }
    X=nxt+ghostx*2;
    Y=nyt+ghosty*2;
    Z=nzt+ghostz*2;
    /*
    for(int ii=0; ii<nxt+ghostx*2; ii++)
    for(int jj=0; jj<nyt+ghosty*2; jj++)
    for(int kk=0; kk<nzt+ghostz*2; kk++)
    {
       u1[ii][jj][kk] = -99;
    }*/

    cleanChunk((float*)u1);

    //for debug info output
    char info[1000000]="\0", tmp[200]="\0";

    //Multi Write test
    time(&start);
    t1=MPI_Wtime();

    for(int i=0; i< TIMES; i++)
    {
        sprintf(tmp, "%d) In main(), %d time write test.\n",rank, i);
        strcpy(info,tmp);
        //MPI_Barrier(Comm);
        unsigned int num =0;
        for(int kk=ghostz; kk<nzt+ghostz*2-ghostz; kk++)
        {
            sum =0;
            for(int jj=ghosty; jj<nyt+ghosty*2-ghosty; jj++)
            for(int ii=ghostx; ii<nxt+ghostx*2-ghostx; ii++)
            {
               if(size>=16)
               {
                    //u1[ii][jj][kk] = rank;
                    //unsigned long p = kk *(X*Y)+ jj*X + ii;
                    unsigned long p = ii *(Y*Z)+ jj*Z + kk;
                    u1[p] = (rank+1) * 100 + i;
               }else
               {   //u1[ii][jj][kk] =(float)(i)*(nx*ny*nz)+ (float)(kk-ghostz)*(nx*ny)+ (float)(rank)*(nxt*nyt)+ (float)(sum++);

                   num = coords[2]*(nx*ny*nzt); //Z director
                   num += (kk-ghostz) *(nx*ny);
                   num += coords[0]*(nx*nyt);
                   num += coords[1]*nxt;
                   num += (jj-ghosty)*nx ;
                   num += (ii-ghostx);
                   num += i*(nx*ny*nz);
                   //u1[ii][jj][kk] = num;
                   //unsigned long p = kk *(X*Y)+ jj*X + ii;
                   unsigned long p = ii *(Y*Z)+ jj*Z + kk;
                   u1[p] = num;
               }
            }
        }

        sprintf(tmp,"%d )Grid:\n",rank);
        sprintf(tmp,"%d) ",rank);
        strcpy(info,tmp);

        for(int kk=ghostz; kk<nzt+ghostz; kk++)
        for(int jj=ghosty; jj<nyt+ghosty; jj++)
        {
           sprintf(tmp,"%3d )%3dx%3d=",rank, kk,jj);
           //strcat(info,tmp);

           for(int ii=ghostx; ii<nxt+ghostx; ii++)
           {
               unsigned long p = ii *(Y*Z)+ jj*Z + kk;
               sprintf(tmp," %4.1f ",u1[p]);
               //sprintf(tmp," %4.1f ",u1[ii][jj][kk]);
              //strcat(info,tmp);
           }
           sprintf(tmp,"\n");
           //strcat(info,tmp);
        }

        //seismIO Write vlm test...
        seism_write(&seism_filex2, u1, &err);  //write volume

        //seismIO Write sfc test...
        seism_write(&seism_filex, u1, &err); //write surface


        //if(rank==0) printf("===============================\n");


    }

    t2=MPI_Wtime();
    MPI_Barrier(MC1);
    time(&end);
    t1 = t2-t1;
    MPI_Allreduce( &t1, &dtime, 1, MPI_DOUBLE, MPI_MAX, MC1);

    seism_file_close(&seism_filex, &err);
    seism_file_close(&seism_filex2, &err);

    sprintf(tmp, "%d) rank, Write Times: %ld(c native) %.1f(MPI_MAX) seconds.\n",rankmc1, end-start, dtime);
    strcpy(info,tmp);
    if(rank%4000==0) printf("%s", info);

    //Read test
    if(nxt<0)
    {
        fmode[0]='r';

        //read surface file
        seism_file_open(sxgro, fmode, &write_step, fdata, &seism_regGridID4_read, &seism_filex4_read, &err);
        MPI_Barrier(Comm);
        if(err !=0)
        {
            printf("SEISM ERROR! Open %s failed.", sxgro);
            return 1;
        }
        else{
            if(rankmc1 ==0) printf("%d) rank, %s opened for reading surface.\n", rankmc1, sxgro);
        }

        cleanChunk((float*)u1);

        //seism_readMesh(sxgro, &nvar, u1, &order, &err);

        seism_read(&seism_filex4_read, u1, &err);
        MPI_Barrier(Comm);
        printChunk(u1, nxt, nyt, nzt, ghostx, ghosty, ghostz, rank);

        seism_file_close(&seism_filex4_read, &err);


        //read volumn file
        seism_file_open(sxgro2, fmode, &write_step, fdata, &seism_regGridID3_read, &seism_filex3_read, &err);
        MPI_Barrier(Comm);
        if(err !=0)
        {
            printf("SEISM ERROR! Open %s failed.", sxgro2);
            return 1;
        }
        else{
            if(rankmc1 ==0) printf("%d) rank, %s opened for reading volumn.\n", rankmc1, sxgro2);
        }

        cleanChunk((float*)u1);


        seism_read(&seism_filex3_read, u1, &err);
        MPI_Barrier(Comm);
        printChunk(u1, nxt, nyt, nzt, ghostx, ghosty, ghostz, rank);

        seism_file_close(&seism_filex3_read, &err);
    }

    if(rankmc1 == 0) printf("%d) write and read done.\n",rank);
    if(u1!=NULL)
    {
      free(u1);
      u1 = NULL;
    }
    seism_finalize(&err);
    MPI_Barrier(Comm);
}

int main(int argc, char* argv[])
{

    MPI_Comm Comm,MC1;

    int rank, size, err;
    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &Comm);
    MPI_Comm_rank(Comm, &rank);
    MPI_Comm_size(Comm, &size);

    int cores = 1728, chunksize = 2;
    if(cores != size) cores = size;
    chunksize = 2;
    if(cores >=24) chunksize = 8;
    //chunksize = 112;
    test(Comm, rank, size, 0, cores, chunksize);
    test(Comm, rank, size, 1, cores, chunksize);
    test(Comm, rank, size, 2, cores, chunksize);
    test(Comm, rank, size, 3, cores, chunksize);
    MPI_Finalize();
    return 0;
}
