#ifndef CL_KERNEL_SOURCE
#define CL_KERNEL_SOURCE
#include<map>
std::string kernel = "\n\
// required definitions\n\
//#define DATATYPE int2\n\
//#define STRIDE_ORDER 18\n\
//#define GRANULARITY_ORDER 4\n\
//#define OPT_ZEROCOPY not needed\n\
//#define OPT_NOTDUMMY\n\
\n\
#define CMD_HELPER(FUNC, NAME) FUNC ## _ ## NAME\n\
#define CMD(FUNC, NAME) CMD_HELPER(FUNC, NAME)\n\
\n\
#define STRIDE (1 << STRIDE_ORDER)\n\
#define GRANULARITY (1 << GRANULARITY_ORDER) \n\
\n\
int reduce_int(int v)    { return v; }\n\
int reduce_int2(int2 v)  { return v.x+v.y; }\n\
int reduce_int4(int4 v)  { return v.x+v.y+v.z+v.w; }\n\
int reduce_int8(int8 v)  { return v.s0+v.s1+v.s2+v.s3+v.s4+v.s5+v.s6+v.s7; }\n\
int reduce_int16(int16 v){ return v.s0+v.s1+v.s2+v.s3+v.s4+v.s5+v.s6+v.s7+v.s8+v.s9+v.sA+v.sB+v.sC+v.sD+v.sE+v.sF; }\n\
\n\
__kernel void initialize(__global DATATYPE *data, const uint n, const int v) {\n\
	const uint id = get_global_id(0);\n\
	const uint stride = get_global_size(0);\n\
	uint offset = 0;\n\
	while( id+offset<n ){\n\
		data[id+offset] = (DATATYPE)v;\n\
		offset += stride;\n\
	}\n\
}\n\
\n\
__kernel void kernel1(__global DATATYPE *data, const uint n) {\n\
	// Get our global thread ID\n\
	const uint id = get_global_id(0);\n\
	const uint  low_order_id = id & ( (uint)(STRIDE - 1));\n\
	const uint high_order_id = id & (~(uint)(STRIDE - 1));\n\
	const uint index = (high_order_id << GRANULARITY_ORDER) | low_order_id;\n\
	const int localid = get_local_id(0);\n\
	const int group_size = get_local_size(0);\n\
/*	__local uint lcount;\n\
	barrier(CLK_LOCAL_MEM_FENCE);\n\
	if( localid==0 )\n\
		lcount = 0;*/\n\
	// Make sure we do not go out of bounds\n\
//	DATATYPE tmp = (DATATYPE)0;\n\
//	#pragma unroll\n\
	for(int i=0; i<GRANULARITY; i++){\n\
		DATATYPE v = data[index+i*STRIDE];\n\
		int v_sum = CMD(reduce, DATATYPE)(v);\n\
		if( v_sum>0 )\n\
			data[index+i*STRIDE] = (DATATYPE)0;\n\
		//tmp = tmp+data[index+i*STRIDE];\n\
		//data[index+i*STRIDE] = (DATATYPE)(index+i*STRIDE);\n\
	}\n\
/*	if( count )\n\
		atomic_add(&lcount, count);\n\
	barrier(CLK_LOCAL_MEM_FENCE);\n\
	// Atomic reduce in global memory\n\
	if(localid==0 && lcount){\n\
		atomic_add((__global int*)result, lcount);\n\
	}*/\n\
}\n\
"; 
std::map<std::string,std::string> source_map = {
std::make_pair("kernel.cl",kernel),
};
#endif