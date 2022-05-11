MSTRINGIFY(
\n
\n__kernel void mm2_kernel1(__global float *A, __global float *B, __global float *C, int ni, int nj, int nk) 
\n{    
\n	int j = get_global_id(0);
\n	int i = get_global_id(1);
\n	
\n	if ((i < ni) && (j < nj))
\n	{ 
\n		int k;
\n		for (k = 0; k < nk; k++)
\n		{
\n			C[i * nj + j] += A[i * nk + k] * B[k * nj + j];
\n		}
\n	}
\n}
\n
\n__kernel void mm2_kernel2(__global float *C, __global float *D, __global float *E, int ni, int nl, int nj) 
\n{    
\n	int j = get_global_id(0);
\n	int i = get_global_id(1);
\n	
\n	if ((i < ni) && (j < nl))
\n	{ 
\n		int k;
\n		for (k = 0; k < nj; k++)
\n		{
\n			E[i * nl + j] += C[i * nj + k] * D[k * nl + j];
\n		}
\n	}
\n} 
);
