MSTRINGIFY(
\n/**
\n * gemm.cl: This file is part of the PolyBench/GPU 1.0 test suite.
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
\n__kernel void gemm(__global DATA_TYPE *a, __global DATA_TYPE *b, __global DATA_TYPE *c, DATA_TYPE alpha, DATA_TYPE beta, int ni, int nj, int nk) 
\n{
\n    	int j = get_global_id(0);
\n	int i = get_global_id(1);
\n
\n	if ((i < ni) && (j < nj))
\n	{	
\n		c[i * nj + j] *= beta;
\n		int k;
\n		for(k=0; k < nk; k++)
\n		{
\n			c[i * nj + j] += alpha * a[i * nk + k] * b[k * nj +j];
\n		}
\n	}
\n}
);
