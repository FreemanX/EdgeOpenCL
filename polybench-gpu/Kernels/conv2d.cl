MSTRINGIFY(
\n/**
\n * 2DConvolution.cl: This file is part of the PolyBench/GPU 1.0 test suite.
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
\n/* Can switch DATA_TYPE between float and double */
\ntypedef float DATA_TYPE;
\n
\n__kernel void Convolution2D_kernel(__global DATA_TYPE *A, __global DATA_TYPE *B, int ni, int nj) 
\n{
\n	int j = get_global_id(0);
\n	int i = get_global_id(1);
\n
\n	DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;
\n	c11 = +0.2;  c21 = +0.5;  c31 = -0.8;
\n	c12 = -0.3;  c22 = +0.6;  c32 = -0.9;
\n	c13 = +0.4;  c23 = +0.7;  c33 = +0.10;
\n	if ((i < (ni-1)) && (j < (nj - 1)) && (i > 0) && (j > 0))
\n	{
\n		B[i*nj + j] =  c11 * A[(i - 1) * nj + (j - 1)]  + c21 * A[(i - 1) * nj + (j + 0)] + c31 * A[(i - 1) * nj + (j + 1)] 
\n		      + c12 * A[(i + 0) * nj + (j - 1)]  + c22 * A[(i + 0) * nj + (j + 0)] + c32 * A[(i + 0) * nj + (j + 1)]
\n		      + c13 * A[(i + 1) * nj + (j - 1)]  + c23 * A[(i + 1) * nj + (j + 0)] + c33 * A[(i + 1) * nj + (j + 1)];
\n	}
\n}
);
