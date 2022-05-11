MSTRINGIFY(
\n/**
\n * covariance.cl: This file is part of the PolyBench/GPU 1.0 test suite.
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
\n
\n__kernel void mean_kernel(__global DATA_TYPE *mean, __global DATA_TYPE *data, DATA_TYPE float_n, int m, int n) 
\n{    	 
\n	int j = get_global_id(0) + 1;
\n	
\n	if ((j >= 1) && (j < (m+1)))
\n	{
\n		mean[j] = 0.0;
\n
\n		int i;
\n		for(i = 1; i < (n+1); i++)
\n		{
\n			mean[j] += data[i * (m+1) + j];
\n		}
\n		mean[j] /= (DATA_TYPE)float_n;
\n	}
\n}
\n
\n__kernel void reduce_kernel(__global DATA_TYPE *mean, __global DATA_TYPE *data, int m, int n) 
\n{
\n	int j = get_global_id(0) + 1;    
\n	int i = get_global_id(1) + 1;
\n
\n	if ((i >= 1) && (i < (n+1)) && (j >= 1) && (j < (m+1)))
\n	{
\n		data[i * (m+1) + j] -= mean[j];	
\n	}
\n}
\n
\n
\n__kernel void covar_kernel(__global DATA_TYPE *symmat, __global DATA_TYPE *data, int m, int n) 
\n{
\n	int j1 = get_global_id(0) + 1;
\n	int i, j2;
\n
\n	if ((j1 >= 1) && (j1 < (m+1)))
\n	{
\n		for (j2 = j1; j2 < (m+1); j2++)
\n		{		
\n	      		symmat[j1*(m+1) + j2] = 0.0;
\n			for(i = 1; i < (n+1); i++)
\n			{
\n				symmat[j1 * (m+1) + j2] += data[i *(m+1) + j1] * data[i *(m+1) + j2];
\n			}
\n			symmat[j2 * (m+1) + j1] = symmat[j1 * (m+1) + j2];
\n		}
\n	}
\n} 
);
