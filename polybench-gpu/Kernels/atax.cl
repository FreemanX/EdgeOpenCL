MSTRINGIFY(
\n/**
\n * atax.cl: This file is part of the PolyBench/GPU 1.0 test suite.
\n *
\n *
\n * Contact: Scott Grauer-Gray <sgrauerg@gmail.com>
\n * Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
\n * Web address: http://www.cse.ohio-state.edu/~pouchet/software/polybench/GPU
\n */
\n
\n#if defined(cl_khr_fp64)  // Khronos extension available?
\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable
\n#elif defined(cl_amd_fp64)  // AMD extension available?
\n#pragma OPENCL EXTENSION cl_amd_fp64 : enable
\n#endif
\n
\ntypedef float DATA_TYPE;
\n
\n__kernel void atax_kernel1(__global DATA_TYPE *A, __global DATA_TYPE *x, __global DATA_TYPE *tmp, int nx, int ny) {
\n    
\n	int i = get_global_id(0);
\n
\n	if (i < nx)
\n	{
\n		int j;
\n		for(j=0; j < ny; j++)
\n		{
\n			tmp[i] += A[i * ny + j] * x[j];
\n		}
\n	}
\n}
\n
\n__kernel void atax_kernel2(__global DATA_TYPE *A, __global DATA_TYPE *y, __global DATA_TYPE *tmp, int nx, int ny) {
\n    
\n	int j = get_global_id(0);
\n
\n	if (j < ny)
\n	{
\n		int i;
\n		for(i=0; i < nx; i++)
\n		{
\n			y[j] += A[i * ny + j] * tmp[i];
\n		}
\n	}
\n} 
);
