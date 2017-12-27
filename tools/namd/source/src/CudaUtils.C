
#include <stdio.h>
#include "common.h"
#include "charm++.h"
#include "CudaUtils.h"

#ifdef NAMD_CUDA

void cudaDie(const char *msg, cudaError_t err) {
  char host[128];
#ifdef NOHOSTNAME
  sprintf(host,"physical node %d", CmiPhysicalNodeID(CkMyPe()));
#else
  gethostname(host, 128);  host[127] = 0;
#endif
  char devstr[128] = "";
  int devnum;
  if ( cudaGetDevice(&devnum) == cudaSuccess ) {
    sprintf(devstr, " device %d", devnum);
  }
  char errmsg[1024];
  if (err == cudaSuccess) {
    sprintf(errmsg,"CUDA error %s on Pe %d (%s%s)", msg, CkMyPe(), host, devstr);
  } else {
    sprintf(errmsg,"CUDA error %s on Pe %d (%s%s): %s", msg, CkMyPe(), host, devstr, cudaGetErrorString(err));    
  }
  NAMD_die(errmsg);
}

//----------------------------------------------------------------------------------------
//
// Deallocate gpu memory
// pp = memory pointer
//
void deallocate_device_T(void **pp) {
  
  if (*pp != NULL) {
    cudaCheck(cudaFree((void *)(*pp)));
    *pp = NULL;
  }

}

//----------------------------------------------------------------------------------------
//
// Deallocate page-locked host memory
// pp = memory pointer
//
void deallocate_host_T(void **pp) {
  
  if (*pp != NULL) {
    cudaCheck(cudaFreeHost((void *)(*pp)));
    *pp = NULL;
  }

}

//----------------------------------------------------------------------------------------
//
// Allocate & re-allocate device memory
// pp = memory pointer
// curlen = current length of the array
// newlen = new required length of the array
// fac = extra space allocation factor: in case of re-allocation new length will be fac*newlen
//
// returns true if reallocation happened
//
bool reallocate_device_T(void **pp, int *curlen, const int newlen, const float fac, const size_t sizeofT) {

  if (*pp != NULL && *curlen < newlen) {
    cudaCheck(cudaFree((void *)(*pp)));
    *pp = NULL;
  }

  if (*pp == NULL) {
    if (fac > 1.0f) {
      *curlen = (int)(((double)(newlen))*(double)fac);
    } else {
      *curlen = newlen;
    }
    cudaCheck(cudaMalloc(pp, sizeofT*(*curlen)));
    return true;
  }

  return false;
}

//----------------------------------------------------------------------------------------
//
// Allocate gpu memory
// pp = memory pointer
// len = length of the array
//
void allocate_device_T(void **pp, const int len, const size_t sizeofT) {
  cudaCheck(cudaMalloc(pp, sizeofT*len));
}

void copy_DtoD_T(const void *d_src, void *d_dst, const int array_len, const size_t sizeofT) {
  cudaCheck(cudaMemcpy(d_dst, d_src, sizeofT*array_len, cudaMemcpyDeviceToDevice));
}

//----------------------------------------------------------------------------------------
//
// Allocate & re-allocate page-locked host memory, preserves content
//
bool resize_device_T(void **pp, int *curlen, const int cur_size, const int new_size,
        const float fac, const size_t sizeofT) {

  void *old = NULL;  

  if (*pp != NULL && *curlen < new_size) {
    allocate_device_T(&old, cur_size, sizeofT);
    copy_DtoD_T(*pp, old, cur_size, sizeofT);
    cudaCheck(cudaDeviceSynchronize());       //Make sure D-D copy is done
    cudaCheck(cudaFree((void *)(*pp)));
    *pp = NULL;
  }

  if (*pp == NULL) {
    if (fac > 1.0f) {
      *curlen = (int)(((double)(new_size))*(double)fac);
    } else {
      *curlen = new_size;
    }
    allocate_device_T(pp, *curlen, sizeofT);
    if (old != NULL) {
      copy_DtoD_T(old, *pp, cur_size, sizeofT);
      cudaCheck(cudaDeviceSynchronize());       //Make sure D-D copy is done
      deallocate_device_T(&old);
    }
    return true;
  }

  return false;
}

//----------------------------------------------------------------------------------------
//
// Allocate & re-allocate page-locked host memory
// pp = memory pointer
// curlen = current length of the array
// newlen = new required length of the array
// fac = extra space allocation factor: in case of re-allocation new length will be fac*newlen
// flag = allocation type:
//        cudaHostAllocDefault = default type, emulates cudaMallocHost
//        cudaHostAllocMapped  = maps allocation into CUDA address space
//
// returns true if reallocation happened
//
bool reallocate_host_T(void **pp, int *curlen, const int newlen, 
		       const float fac, const unsigned int flag, const size_t sizeofT) {

  if (*pp != NULL && *curlen < newlen) {
    cudaCheck(cudaFreeHost((void *)(*pp)));
    *pp = NULL;
  }

  if (*pp == NULL) {
    if (fac > 1.0f) {
      *curlen = (int)(((double)(newlen))*(double)fac);
    } else {
      *curlen = newlen;
    }
    cudaCheck(cudaHostAlloc(pp, sizeofT*(*curlen), flag));
    return true;
  }

  return false;
}

//----------------------------------------------------------------------------------------
//
// Allocate & re-allocate page-locked host memory, preserves content
//
bool resize_host_T(void **pp, int *curlen, const int cur_size, const int new_size,
       const float fac, const unsigned int flag, const size_t sizeofT) {

  char *old = NULL;

  if (*pp != NULL && *curlen < new_size) {
    old = new char[cur_size*sizeofT];
    memcpy(old, *pp, cur_size*sizeofT);
    cudaCheck(cudaFreeHost((void *)(*pp)));
    *pp = NULL;
  }

  if (*pp == NULL) {
    if (fac > 1.0f) {
      *curlen = (int)(((double)(new_size))*(double)fac);
    } else {
      *curlen = new_size;
    }
    cudaCheck(cudaHostAlloc(pp, sizeofT*(*curlen), flag));
    if (old != NULL) {
      memcpy(*pp, old, cur_size*sizeofT);
      delete [] old;
    }
    return true;
  }

  return false;
}

#endif // NAMD_CUDA
