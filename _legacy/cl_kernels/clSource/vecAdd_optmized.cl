__kernel void matAdd(__global const float* A,
	__global const float* B, 
	__global float* C,
	const int N) {
	int id = get_global_id(0);
	int idx = mul24(id, 4);
	float4 tmpA = vload4(idx, A);
	float4 tmpB = vload4(idx, B);
	(*((__global float4*)&C[idx])) = (tmpA+tmpB);
}

__kernel void matAdd2(
	__global float* matrixA,
	__global float* matrixB,
	__global float* MatrixSum,
	const int rows,
	const int cols)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	int offset = mul24(j, cols);
	int index = mad24(i, 4, offset);

	float4 tmpA = (*((__global float4*)&matrixA[index])); 
	vload and vstore can be used in here
	float4 tmpB = (*((__global float4*)&matrixB[index]));
	(*((__global float4*)&MatrixSum[index])) = (tmpA+tmpB);

}