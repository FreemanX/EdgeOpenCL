MSTRINGIFY(
\n/**
\n * gesummv.cl: This file is part of the PolyBench/GPU 1.0 test suite.
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
\n	
\n__kernel void gesummv_kernel(__global DATA_TYPE *a, __global DATA_TYPE *b, __global DATA_TYPE *x, __global DATA_TYPE *y, __global DATA_TYPE *tmp, DATA_TYPE alpha, DATA_TYPE beta, int n) 
\n{    
\n	int i = get_global_id(0);
\n
\n	if (i < n)
\n	{
\n		int j;
\n		for(j = 0; j < n; j++)
\n		{	
\n			tmp[i] += a[i * n + j] * x[j];
\n			y[i] += b[i * n + j] * x[j];
\n		}
\n		y[i] = alpha * tmp[i] + beta * y[i];
\n	}
\n}
);
