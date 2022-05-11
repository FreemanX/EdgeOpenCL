MSTRINGIFY(
\n/**
\n * gramschmidt.cl: This file is part of the PolyBench/GPU 1.0 test suite.
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
\n__kernel void gramschmidt_kernel1(__global DATA_TYPE *a, __global DATA_TYPE *r, __global DATA_TYPE *q, int k, int m, int n)
\n{
\n	int tid = get_global_id(0);
\n	
\n	if (tid == 0)
\n	{
\n		DATA_TYPE nrm = 0.0;
\n		int i;
\n		for (i = 0; i < m; i++)
\n		{
\n			nrm += a[i * n + k] * a[i * n + k];
\n		}
\n      		r[k * n + k] = sqrt(nrm);
\n	}
\n}
\n
\n
\n__kernel void gramschmidt_kernel2(__global DATA_TYPE *a, __global DATA_TYPE *r, __global DATA_TYPE *q, int k, int m, int n)
\n{
\n	int i = get_global_id(0);
\n
\n        if (i < m)
\n	{	
\n		q[i * n + k] = a[i * n + k] / r[k * n + k];
\n	}
\n}
\n
\n
\n__kernel void gramschmidt_kernel3(__global DATA_TYPE *a, __global DATA_TYPE *r, __global DATA_TYPE *q, int k, int m, int n)
\n{
\n	int j = get_global_id(0);
\n
\n	if ((j > k) && (j < n))
\n	{
\n		r[k*n + j] = 0.0;
\n
\n		int i;
\n		for (i = 0; i < m; i++)
\n		{
\n			r[k*n + j] += q[i*n + k] * a[i*n + j];
\n		}
\n		
\n		for (i = 0; i < m; i++)
\n		{
\n			a[i*n + j] -= q[i*n + k] * r[k*n + j];
\n		}
\n	}
\n}
);
