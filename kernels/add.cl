__kernel void vecAdd(
__global float* A,
__global float* B,
__global float* C,
const unsigned int N) {
	int i = get_global_id(0);
	if (i < N) {
		C[i] = A[i] + B[i];
	}
}

__kernel void matAdd(
__global float* matrixA,
__global float* matrixB,
__global float* MatrixSum,
const int rows,
const int cols) {
	int i = get_global_id(0);
	int j = get_global_id(1);
	int offset = mul24(j, cols);
	int index = mad24(i, 4, offset);

	float4 tmpA = (*((__global float4*)&matrixA[index]));
	float4 tmpB = (*((__global float4*)&matrixB[index]));
	(*((__global float4*)&MatrixSum[index])) += (tmpA+tmpB);
}
