MSTRINGIFY(
\n
\n __kernel void vecAdd(
\n __global float* A,
\n __global float* B,
\n __global float* C,
\n const unsigned int N) {
\n 	int i = get_global_id(0);
\n 	if (i < N) {
\n 		C[i] = A[i] + B[i];
\n 	}
\n }
\n
\n __kernel void matAdd(
\n __global float* matrixA,
\n __global float* matrixB,
\n __global float* MatrixSum,
\n const int rows,
\n const int cols) {
\n 	int i = get_global_id(0);
\n 	int j = get_global_id(1);
\n 	int offset = mul24(j, cols);
\n 	int index = mad24(i, 4, offset);
\n
\n 	float4 tmpA = (*((__global float4*)&matrixA[index]));
\n 	float4 tmpB = (*((__global float4*)&matrixB[index]));
\n 	(*((__global float4*)&MatrixSum[index])) += (tmpA+tmpB);
\n }
\n
);
