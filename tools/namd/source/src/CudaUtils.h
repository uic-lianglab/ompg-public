#ifndef CUDAUTILS_H
#define CUDAUTILS_H

#ifdef NAMD_CUDA
#include <cuda.h>
#include <cuda_runtime.h>

//
// Cuda static assert, copied from Facebook FFT sources. Remove once nvcc has c++11
//
template <bool>
struct CudaStaticAssert;

template <>
struct CudaStaticAssert<true> {
};

#define cuda_static_assert(expr) \
  (CudaStaticAssert<(expr) != 0>())

void cudaDie(const char *msg, cudaError_t err=cudaSuccess);

//
// Error checking wrapper for CUDA
//
#define cudaCheck(stmt) do {                                 \
	cudaError_t err = stmt;                            \
  if (err != cudaSuccess) {                          \
  	char msg[128];	\
	  sprintf(msg, "%s in file %s, function %s\n", #stmt,__FILE__,__FUNCTION__); \
	  cudaDie(msg, err); \
  }                                                  \
} while(0)

void deallocate_device_T(void **pp);
//----------------------------------------------------------------------------------------
//
// Deallocate gpu memory
// pp = memory pointer
//
template <class T>
void deallocate_device(T **pp) {
  deallocate_device_T((void **)pp);
}
//----------------------------------------------------------------------------------------

bool resize_device_T(void **pp, int *curlen, const int cur_size, const int new_size,
        const float fac, const size_t sizeofT);

//----------------------------------------------------------------------------------------
//
// Allocate & re-allocate GPU memory, preserves content
// pp = memory pointer
// curlen = current length of the array
// cur_size = current size of the data content
// new_size = new size
// fac = extra space allocation factor: in case of re-allocation new length will be fac*new_size
//
template <class T>
bool resize_device(T **pp, int *curlen, const int cur_size, const int new_size, const float fac) {
  return resize_device_T((void **)pp, curlen, cur_size, new_size, fac, sizeof(T));
}

bool reallocate_device_T(void **pp, int *curlen, const int newlen, const float fac, const size_t sizeofT);
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
template <class T>
bool reallocate_device(T **pp, int *curlen, const int newlen, const float fac=1.0f) {
  return reallocate_device_T((void **)pp, curlen, newlen, fac, sizeof(T));
}
//----------------------------------------------------------------------------------------
bool reallocate_host_T(void **pp, int *curlen, const int newlen, const float fac, 
		       const unsigned int flag, const size_t sizeofT);
//----------------------------------------------------------------------------------------
//
// Allocate & re-allocate pinned host memory
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
template <class T>
bool reallocate_host(T **pp, int *curlen, const int newlen,
		     const float fac=1.0f, const unsigned int flag=cudaHostAllocDefault) {
  return reallocate_host_T((void **)pp, curlen, newlen, fac, flag, sizeof(T));
}

bool resize_host_T(void **pp, int *curlen, const int cur_size, const int new_size,
       const float fac, const unsigned int flag, const size_t sizeofT);

//----------------------------------------------------------------------------------------
//
// Allocate & re-allocate host memory, preserves content
// pp = memory pointer
// curlen = current length of the array
// cur_size = current size of the data content
// new_size = new size
// fac = extra space allocation factor: in case of re-allocation new length will be fac*new_size
// flag = allocation type:
//        cudaHostAllocDefault = default type, emulates cudaMallocHost
//        cudaHostAllocMapped  = maps allocation into CUDA address space
//
template <class T>
bool resize_host(T **pp, int *curlen, const int cur_size, const int new_size,
  const float fac=1.0f, const unsigned int flag=cudaHostAllocDefault) {
  return resize_host_T((void **)pp, curlen, cur_size, new_size, fac, flag, sizeof(T));
}

void deallocate_host_T(void **pp);
//----------------------------------------------------------------------------------------
//
// Deallocate page-locked host memory
// pp = memory pointer
//
template <class T>
void deallocate_host(T **pp) {
  deallocate_host_T((void **)pp);
}
//----------------------------------------------------------------------------------------

#endif // NAMD_CUDA

#endif // CUDAUTILS_H