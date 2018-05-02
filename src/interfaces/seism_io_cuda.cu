/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
/*************************************************************************
	> File Name: seism_io_cuda.cu
	> Author: huizhou
	> Mail: wbzhui@gmail.com
	> Created Time: Wed 07 Sep 2016 10:47:41 AM PDT
 ************************************************************************/

#include<stdio.h>
#include<cuda.h>
#include<cuda_runtime.h>
#include<assert.h>
#include "seism_io_globals.h"
#include "seism_io_consts.h"

int bufferfromGPU = 0;

/*
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    }
}
*/

void *convertCudaPointer(FileProfile *fp, void *buffer, int *err)
{

    cudaError_t code;

    //cudaDeviceProp prop;
    //code = cudaGetDeviceProperties(&prop , 0);

    cudaPointerAttributes attributes;
    code =cudaPointerGetAttributes (&attributes, buffer);
    if (code != cudaSuccess)
    {
        fprintf(stderr,"cudaPointerGetAttributes error: %s.\n", cudaGetErrorString(code));
        *err = -1;
        return NULL;
    }

    //if buffer comes from GPU Space
    if(attributes.memoryType==cudaMemoryTypeDevice)
    {
        bufferfromGPU = 1;
        PointSet *ps = fp->pointSet;
        if(!ps)
        {
            printf("ERROR! SEISM_IO: no PS in FP!\n");
            *err = -2;
            return NULL;
        }

    	//void* d_data = (void*)cudaMallocHost(sizeof(fp->datatype)*(ps->nprec)*(fp->write_step));
        void* m_data = (void*)malloc(sizeof(fp->datatype)*(ps->nprec)*(fp->write_step));
        //if (code != cudaSuccess)
        if (m_data != NULL)
        {
            //fprintf(stderr,"cudaMallocHost error: %s.\n", cudaGetErrorString(code));
            fprintf(stderr,"Malloc error.\n");
            *err = -3;
            return NULL;
        }

        code = cudaMemcpy(m_data, buffer, sizeof(fp->datatype)*(ps->nprec)*(fp->write_step), cudaMemcpyDeviceToHost);
        if (code != cudaSuccess)
        {
            fprintf(stderr,"cudaMemcpy error: %s.\n", cudaGetErrorString(code));
            //cudaFree(buffer);
            free(buffer);
            buffer = NULL;
            *err = -4;
            return NULL;
        }
        *err = 0;
        return m_data;

    }
    //buffer comes from Host
    else
    {
        bufferfromGPU = 0;
        return buffer;
    }

}

void freeCudaPointer(void *buffer)
{
    if(bufferfromGPU = 1 && buffer != NULL)
    {
        //cudaFree(buffer);
        free(buffer);
        buffer = NULL;
    }
}
