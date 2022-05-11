MSTRINGIFY(
\n/**
\n * syr.cl: This file is part of the PolyBench/GPU 1.0 test suite.
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
\n__kernel void syr2k_kernel(__global DATA_TYPE *a, __global DATA_TYPE *b, __global DATA_TYPE *c, DATA_TYPE alpha, DATA_TYPE beta, int m, int n) 
\n{    
\n   	int j = get_global_id(0);
\n	int i = get_global_id(1);
\n
\n	if ((i < n) && (j < n))
\n	{
\n		c[i * n + j] *= beta;
\n		
\n		int k;
\n		for(k = 0; k < m; k++)
\n		{
\n			c[i * n + j] += alpha * a[i * m + k] * b[j * m + k] + alpha * b[i * m + k] * a[j * m + k];
\n		}
\n	}
\n}
);
