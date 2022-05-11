//
// Created by pfxu on 7/19/19.
//

#include "benchmark.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

static const char *PROGRAM_SOURCE[] = {
        "#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable\n",
        "\n",
        "__kernel void histo_prescan_kernel (__global unsigned int* input, int size, __global unsigned int* minmax)\n",
        "{\n",
        "\n",
        "    __local float Avg[PRESCAN_THREADS];\n",
        "    __local float StdDev[PRESCAN_THREADS];\n",
        "\n",
        "    int threadIdxx = get_local_id(0);\n",
        "    int blockDimx = get_local_size(0);\n",
        "    int blockIdxx = get_group_id(0);\n",
        "    int stride = size/(get_num_groups(0));\n",
        "    int addr = blockIdxx*stride+threadIdxx;\n",
        "    int end = blockIdxx*stride + stride/8; // Only sample 1/8th of the input data\n",
        "\n",
        "    float avg = 0.0;\n",
        "    unsigned int count = 0;\n",
        "    while (addr < end){\n",
        "      avg += input[addr];\n",
        "      count++;\n",
        "	  addr += blockDimx;\n",
        "    }\n",
        "    avg /= count;\n",
        "    Avg[threadIdxx] = avg;\n",
        "\n",
        "    int addr2 = blockIdxx*stride+threadIdxx;\n",
        "    float stddev = 0;\n",
        "    while (addr2 < end){\n",
        "        stddev += (input[addr2]-avg)*(input[addr2]-avg);\n",
        "        addr2 += blockDimx;\n",
        "    }\n",
        "    stddev /= count;\n",
        "    StdDev[threadIdxx] = sqrt(stddev);\n",
        "\n",
        "#define SUM(stride__) if(threadIdxx < stride__){ Avg[threadIdxx] += Avg[threadIdxx+stride__]; StdDev[threadIdxx] += StdDev[threadIdxx+stride__];}\n",
        "\n",
        "#if (PRESCAN_THREADS >= 32)\n",
        "    for (int stride = PRESCAN_THREADS/2; stride >= 32; stride = stride >> 1){\n",
        "	barrier(CLK_LOCAL_MEM_FENCE);\n",
        "	SUM(stride);\n",
        "    }\n",
        "#endif\n",
        "#if (PRESCAN_THREADS >= 16)\n",
        "    SUM(16);\n",
        "#endif\n",
        "#if (PRESCAN_THREADS >= 8)\n",
        "    SUM(8);\n",
        "#endif\n",
        "#if (PRESCAN_THREADS >= 4)\n",
        "    SUM(4);\n",
        "#endif\n",
        "#if (PRESCAN_THREADS >= 2)\n",
        "    SUM(2);\n",
        "#endif\n",
        "\n",
        "    if (threadIdxx == 0){\n",
        "        float avg = Avg[0]+Avg[1];\n",
        "	avg /= PRESCAN_THREADS;\n",
        "	float stddev = StdDev[0]+StdDev[1];\n",
        "	stddev /= PRESCAN_THREADS;\n",
        "	    atom_min(minmax,((unsigned int)(avg-10*stddev))/(KB*1024));\n",
        "        atom_max(minmax+1,((unsigned int)(avg+10*stddev))/(KB*1024));\n",
        "    }\n",
        "}\n",
        "\n",
        "#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable\n",
        "#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable\n",
        "\n",
        "void testIncrementGlobal (\n",
        "        __global unsigned int *global_histo,\n",
        "        unsigned int sm_range_min,\n",
        "        unsigned int sm_range_max,\n",
        "        const uchar4 sm)\n",
        "{\n",
        "        const unsigned int range = sm.x;\n",
        "        const unsigned int indexhi = sm.y;\n",
        "        const unsigned int indexlo = sm.z;\n",
        "        const unsigned int offset  = sm.w;\n",
        "\n",
        "        if (range < sm_range_min || range > sm_range_max)\n",
        "        {\n",
        "                const unsigned int bin = range * BINS_PER_BLOCK + offset / 8 + (indexlo << 2) + (indexhi << 10);\n",
        "                const unsigned int bin_div2 = bin / 2;\n",
        "                const unsigned int bin_offset = (bin % 2 == 1) ? 16 : 0;\n",
        "\n",
        "                unsigned int old_val = global_histo[bin_div2];\n",
        "                unsigned short old_bin = (old_val >> bin_offset) & 0xFFFF;\n",
        "\n",
        "                if (old_bin < 255)\n",
        "                {\n",
        "                        atom_add (&global_histo[bin_div2], 1 << bin_offset);\n",
        "                }\n",
        "        }\n",
        "}\n",
        "\n",
        "void testIncrementLocal (\n",
        "        __global unsigned int *global_overflow,\n",
        "        __local unsigned int smem[KB][256],\n",
        "        const unsigned int myRange,\n",
        "        const uchar4 sm)\n",
        "{\n",
        "        const unsigned int range = sm.x;\n",
        "        const unsigned int indexhi = sm.y;\n",
        "        const unsigned int indexlo = sm.z;\n",
        "        const unsigned int offset  = sm.w;\n",
        "\n",
        "        if (range == myRange)\n",
        "        {\n",
        "                // Atomically increment shared memory\n",
        "                unsigned int add = (unsigned int)(1 << offset);\n",
        "                unsigned int prev = atom_add (&smem[indexhi][indexlo], add);\n",
        "\n",
        "                // Check if current bin overflowed\n",
        "                unsigned int prev_bin_val = (prev >> offset) & 0x000000FF;\n",
        "\n",
        "                // If there was an overflow, record it and record if it cascaded into other bins\n",
        "                if (prev_bin_val == 0x000000FF)\n",
        "                {\n",
        "                        const unsigned int bin =\n",
        "                                range * BINS_PER_BLOCK +\n",
        "                                offset / 8 + (indexlo << 2) + (indexhi << 10);\n",
        "\n",
        "                        bool can_overflow_to_bin_plus_1 = (offset < 24) ? true : false;\n",
        "                        bool can_overflow_to_bin_plus_2 = (offset < 16) ? true : false;\n",
        "                        bool can_overflow_to_bin_plus_3 = (offset <  8) ? true : false;\n",
        "\n",
        "                        bool overflow_into_bin_plus_1 = false;\n",
        "                        bool overflow_into_bin_plus_2 = false;\n",
        "                        bool overflow_into_bin_plus_3 = false;\n",
        "\n",
        "                        unsigned int prev_bin_plus_1_val = (prev >> (offset +  8)) & 0x000000FF;\n",
        "                        unsigned int prev_bin_plus_2_val = (prev >> (offset + 16)) & 0x000000FF;\n",
        "                        unsigned int prev_bin_plus_3_val = (prev >> (offset + 24)) & 0x000000FF;\n",
        "\n",
        "                        if (can_overflow_to_bin_plus_1 &&        prev_bin_val == 0x000000FF) overflow_into_bin_plus_1 = true;\n",
        "                        if (can_overflow_to_bin_plus_2 && prev_bin_plus_1_val == 0x000000FF) overflow_into_bin_plus_2 = true;\n",
        "                        if (can_overflow_to_bin_plus_3 && prev_bin_plus_2_val == 0x000000FF) overflow_into_bin_plus_3 = true;\n",
        "\n",
        "                        unsigned int bin_plus_1_add;\n",
        "                        unsigned int bin_plus_2_add;\n",
        "                        unsigned int bin_plus_3_add;\n",
        "\n",
        "                        if (overflow_into_bin_plus_1) bin_plus_1_add = (prev_bin_plus_1_val < 0x000000FF) ? 0xFFFFFFFF : 0x000000FF;\n",
        "                        if (overflow_into_bin_plus_2) bin_plus_2_add = (prev_bin_plus_2_val < 0x000000FF) ? 0xFFFFFFFF : 0x000000FF;\n",
        "                        if (overflow_into_bin_plus_3) bin_plus_3_add = (prev_bin_plus_3_val < 0x000000FF) ? 0xFFFFFFFF : 0x000000FF;\n",
        "\n",
        "                                                      atom_add (&global_overflow[bin],  256);\n",
        "                        if (overflow_into_bin_plus_1) atom_add (&global_overflow[bin+1], bin_plus_1_add);\n",
        "                        if (overflow_into_bin_plus_2) atom_add (&global_overflow[bin+2], bin_plus_2_add);\n",
        "                        if (overflow_into_bin_plus_3) atom_add (&global_overflow[bin+3], bin_plus_3_add);\n",
        "                }\n",
        "        }\n",
        "}\n",
        "\n",
        "void clearMemory (__local unsigned int smem[KB][256])\n",
        "{\n",
        "        for (int i = get_local_id(0), blockDimx = get_local_size(0); i < BINS_PER_BLOCK / 4; i += blockDimx)\n",
        "        {\n",
        "                ((__local unsigned int*)smem)[i] = 0;\n",
        "        }\n",
        "}\n",
        "\n",
        "void copyMemory (__global unsigned int *dst, __local unsigned int src[KB][256])\n",
        "{\n",
        "        for (int i = get_local_id(0), blockDimx = get_local_size(0); i < BINS_PER_BLOCK/4; i += blockDimx)\n",
        "        {\n",
        "                atom_add(dst+i*4, (((__local unsigned int*)src)[i] >> 0) & 0xFF );\n",
        "                atom_add(dst+i*4+1, (((__local unsigned int*)src)[i] >> 8) & 0xFF);\n",
        "                atom_add(dst+i*4+2, (((__local unsigned int*)src)[i] >> 16) & 0xFF);\n",
        "                atom_add(dst+i*4+3, (((__local unsigned int*)src)[i] >> 24) & 0xFF);\n",
        "                //dst[i] = ((__local unsigned int*)src)[i];\n",
        "        }\n",
        "}\n",
        "\n",
        "__kernel void histo_main_kernel (\n",
        "        __global uchar4 *sm_mappings,\n",
        "        unsigned int num_elements,\n",
        "        unsigned int sm_range_min,\n",
        "        unsigned int sm_range_max,\n",
        "        unsigned int histo_height,\n",
        "        unsigned int histo_width,\n",
        "        __global unsigned int *global_subhisto,\n",
        "        __global unsigned int *global_histo,\n",
        "        __global unsigned int *global_overflow)\n",
        "{\n",
        "        /* Most optimal solution uses 24 * 1024 bins per threadblock */\n",
        "        __local unsigned int sub_histo[KB][256];\n",
        "\n",
        "        /* Each threadblock contributes to a specific 24KB range of histogram,\n",
        "         * and also scans every N-th line for interesting data.  N = gridDim.x\n",
        "         */\n",
        "        unsigned int blockDimx = get_local_size(0);\n",
        "        unsigned int gridDimx = get_num_groups(0);\n",
        "        unsigned int local_scan_range = sm_range_min + get_group_id(1);\n",
        "        unsigned int local_scan_load = get_group_id(0) * blockDimx + get_local_id(0);\n",
        "\n",
        "        clearMemory (sub_histo);\n",
        "        barrier(CLK_LOCAL_MEM_FENCE); //mem_fence(CLK_GLOBAL_MEM_FENCE);//        __syncthreads();\n",
        "\n",
        "        if (get_group_id(1) == 0 )\n",
        "        {\n",
        "                // Loop through and scan the input\n",
        "                while (local_scan_load < num_elements)\n",
        "                {\n",
        "                        // Read buffer\n",
        "                        uchar4 sm = sm_mappings[local_scan_load];\n",
        "                        local_scan_load += blockDimx * gridDimx;\n",
        "\n",
        "                        // Check input\n",
        "                        testIncrementLocal (\n",
        "                                global_overflow,\n",
        "                                sub_histo,\n",
        "                                local_scan_range,\n",
        "                                sm\n",
        "                        );\n",
        "                        testIncrementGlobal (\n",
        "                                global_histo,\n",
        "                                sm_range_min,\n",
        "                                sm_range_max,\n",
        "                                sm\n",
        "                        );\n",
        "                }\n",
        "        }\n",
        "        else\n",
        "        {\n",
        "                // Loop through and scan the input\n",
        "                while (local_scan_load < num_elements)\n",
        "                {\n",
        "                        // Read buffer\n",
        "                        uchar4 sm = sm_mappings[local_scan_load];\n",
        "                        local_scan_load += blockDimx * gridDimx;\n",
        "\n",
        "                        // Check input\n",
        "                        testIncrementLocal (\n",
        "                                global_overflow,\n",
        "                                sub_histo,\n",
        "                                local_scan_range,\n",
        "                                sm\n",
        "                        );\n",
        "                }\n",
        "        }\n",
        "\n",
        "        // Store sub histogram to global memory\n",
        "        unsigned int store_index = (local_scan_range * BINS_PER_BLOCK);\n",
        "\n",
        "        barrier(CLK_LOCAL_MEM_FENCE); //__syncthreads();\n",
        "        copyMemory (&(global_subhisto[store_index]), sub_histo);\n",
        "}\n",
        "\n",
        "__kernel void histo_final_kernel (\n",
        "    unsigned int sm_range_min,\n",
        "    unsigned int sm_range_max,\n",
        "    unsigned int histo_height,\n",
        "    unsigned int histo_width,\n",
        "    __global unsigned int *global_subhisto,\n",
        "    __global unsigned int *global_histo,\n",
        "    __global unsigned int *global_overflow,\n",
        "    __global unsigned int *final_histo) //final output\n",
        "{\n",
        "    unsigned int blockDimx = get_local_size(0);\n",
        "    unsigned int gridDimx = get_num_groups(0);\n",
        "    unsigned int start_offset = get_local_id(0) + get_group_id(0) * blockDimx;\n",
        "    const ushort4 zero_short  = {0, 0, 0, 0};\n",
        "    const uint4 zero_int      = {0, 0, 0, 0};\n",
        "\n",
        "    unsigned int size_low_histo = sm_range_min * BINS_PER_BLOCK;\n",
        "    unsigned int size_mid_histo = (sm_range_max - sm_range_min +1) * BINS_PER_BLOCK;\n",
        "\n",
        "    /* Clear lower region of global histogram */\n",
        "    for (unsigned int i = start_offset; i < size_low_histo/4; i += gridDimx * blockDimx)\n",
        "    {\n",
        "        ushort4 global_histo_data = ((__global ushort4*)global_histo)[i];\n",
        "        ((__global ushort4*)global_histo)[i] = zero_short;\n",
        "\n",
        "        global_histo_data.x = min (global_histo_data.x, (ushort) 255);\n",
        "        global_histo_data.y = min (global_histo_data.y, (ushort) 255);\n",
        "        global_histo_data.z = min (global_histo_data.z, (ushort) 255);\n",
        "        global_histo_data.w = min (global_histo_data.w, (ushort) 255);\n",
        "\n",
        "        uchar4 final_histo_data = (uchar4) (\n",
        "            (unsigned char) global_histo_data.x,\n",
        "            (unsigned char) global_histo_data.y,\n",
        "            (unsigned char) global_histo_data.z,\n",
        "            (unsigned char) global_histo_data.w\n",
        "        );\n",
        "\n",
        "        ((__global uchar4*)final_histo)[i] = final_histo_data;\n",
        "    }\n",
        "\n",
        "    /* Clear the middle region of the overflow buffer */\n",
        "    for (unsigned int i = (size_low_histo/4) + start_offset; i < (size_low_histo+size_mid_histo)/4; i += gridDimx * blockDimx)\n",
        "    {\n",
        "        uint4 global_histo_data = ((__global uint4*)global_overflow)[i];\n",
        "        ((__global uint4*)global_overflow)[i] = zero_int;\n",
        "\n",
        "        uint4 internal_histo_data = (uint4)(\n",
        "            global_histo_data.x,\n",
        "            global_histo_data.y,\n",
        "            global_histo_data.z,\n",
        "            global_histo_data.w\n",
        "        );\n",
        "\n",
        "        unsigned int bin4in0 = ((__global unsigned int*)global_subhisto)[i*4];\n",
        "        unsigned int bin4in1 = ((__global unsigned int*)global_subhisto)[i*4+1];\n",
        "        unsigned int bin4in2 = ((__global unsigned int*)global_subhisto)[i*4+2];\n",
        "        unsigned int bin4in3 = ((__global unsigned int*)global_subhisto)[i*4+3];\n",
        "\n",
        "        internal_histo_data.x = min (bin4in0, (unsigned int) 255);\n",
        "        internal_histo_data.y = min (bin4in1, (unsigned int) 255);\n",
        "        internal_histo_data.z = min (bin4in2, (unsigned int) 255);\n",
        "        internal_histo_data.w = min (bin4in3, (unsigned int) 255);\n",
        "\n",
        "        uchar4 final_histo_data = (uchar4) (\n",
        "            internal_histo_data.x,\n",
        "            internal_histo_data.y,\n",
        "            internal_histo_data.z,\n",
        "            internal_histo_data.w\n",
        "        );\n",
        "\n",
        "        ((__global uchar4*)final_histo)[i] = final_histo_data;\n",
        "    }\n",
        "\n",
        "    /* Clear the upper region of global histogram */\n",
        "    for (unsigned int i = ((size_low_histo+size_mid_histo)/4) + start_offset; i < (histo_height*histo_width)/4; i += gridDimx * blockDimx)\n",
        "    {\n",
        "        ushort4 global_histo_data = ((__global ushort4*)global_histo)[i];\n",
        "        ((__global ushort4*)global_histo)[i] = zero_short;\n",
        "\n",
        "        global_histo_data.x = min (global_histo_data.x, (ushort) 255);\n",
        "        global_histo_data.y = min (global_histo_data.y, (ushort) 255);\n",
        "        global_histo_data.z = min (global_histo_data.z, (ushort) 255);\n",
        "        global_histo_data.w = min (global_histo_data.w, (ushort) 255);\n",
        "\n",
        "        uchar4 final_histo_data = (uchar4) (\n",
        "            global_histo_data.x,\n",
        "            global_histo_data.y,\n",
        "            global_histo_data.z,\n",
        "            global_histo_data.w\n",
        "        );\n",
        "\n",
        "        ((__global uchar4*)final_histo)[i] = final_histo_data;\n",
        "    }\n",
        "}\n",
        "\n",
        "__kernel void calculateBin (\n",
        "        __const unsigned int bin,\n",
        "        __global uchar4 *sm_mapping)\n",
        "{\n",
        "        unsigned char offset  =  bin        %   4;\n",
        "        unsigned char indexlo = (bin >>  2) % 256;\n",
        "        unsigned char indexhi = (bin >> 10) %  KB;\n",
        "        unsigned char block   =  bin / BINS_PER_BLOCK;\n",
        "\n",
        "        offset *= 8;\n",
        "\n",
        "        uchar4 sm;\n",
        "        sm.x = block;\n",
        "        sm.y = indexhi;\n",
        "        sm.z = indexlo;\n",
        "        sm.w = offset;\n",
        "\n",
        "        *sm_mapping = sm;\n",
        "}\n",
        "\n",
        "__kernel void histo_intermediates_kernel (\n",
        "        __global uint2 *input,\n",
        "        unsigned int height,\n",
        "        unsigned int width,\n",
        "        unsigned int input_pitch,\n",
        "        __global uchar4 *sm_mappings)\n",
        "{\n",
        "        int threadIdxx = get_local_id(0);\n",
        "        int blockDimx = get_local_size(0);\n",
        "        unsigned int line = UNROLL * (get_group_id(0));// 16 is the unroll factor;\n",
        "\n",
        "        __global uint2 *load_bin = input + line * input_pitch + threadIdxx;\n",
        "\n",
        "        unsigned int store = line * width + threadIdxx;\n",
        "        bool skip = (width % 2) && (threadIdxx == (blockDimx - 1));\n",
        "\n",
        "        #pragma unroll\n",
        "        for (int i = 0; i < UNROLL; i++)\n",
        "        {\n",
        "                uint2 bin_value = *load_bin;\n",
        "\n",
        "                calculateBin (\n",
        "                        bin_value.x,\n",
        "                        &sm_mappings[store]\n",
        "                );\n",
        "\n",
        "                if (!skip) calculateBin (\n",
        "                        bin_value.y,\n",
        "                        &sm_mappings[store + blockDimx]\n",
        "                );\n",
        "\n",
        "                load_bin += input_pitch;\n",
        "                store += width;\n",
        "        }\n",
        "}\n",
        "\n",
        "__kernel void histo_intermediates_kernel_compat (\n",
        "        __global uint2 *input,\n",
        "        unsigned int height,\n",
        "        unsigned int width,\n",
        "        unsigned int input_pitch,\n",
        "        __global uchar4 *sm_mappings)\n",
        "{\n",
        "        int threadIdxx = get_local_id(0);\n",
        "        int blockDimx = input_pitch; //get_local_size(0);\n",
        "\n",
        "        int tid2 = get_local_id(0) + get_local_size(0);\n",
        "\n",
        "        unsigned int line = UNROLL * (get_group_id(0));// 16 is the unroll factor;\n",
        "\n",
        "        __global uint2 *load_bin = input + line * input_pitch + threadIdxx;\n",
        "        __global uint2 *load_bin2 = input + line * input_pitch + tid2;\n",
        "\n",
        "        unsigned int store = line * width + threadIdxx;\n",
        "        unsigned int store2 = line * width + tid2;\n",
        "\n",
        "        bool skip = (width % 2) && (threadIdxx == (input_pitch - 1));\n",
        "        bool skip2 = (width % 2) && (tid2 == (input_pitch - 1));\n",
        "\n",
        "        bool does2 = tid2 < input_pitch;\n",
        "\n",
        "        #pragma unroll\n",
        "        for (int i = 0; i < UNROLL; i++)\n",
        "        {\n",
        "                uint2 bin_value = *load_bin;\n",
        "\n",
        "\n",
        "                calculateBin (\n",
        "                        bin_value.x,\n",
        "                        &sm_mappings[store]\n",
        "                );\n",
        "\n",
        "                if (!skip) calculateBin (\n",
        "                        bin_value.y,\n",
        "                        &sm_mappings[store + blockDimx]\n",
        "                );\n",
        "\n",
        "                load_bin += input_pitch;\n",
        "                store += width;\n",
        "\n",
        "                if (does2) {\n",
        "                  uint2 bin_val2 = *load_bin2;\n",
        "\n",
        "                  calculateBin (\n",
        "                        bin_val2.x,\n",
        "                        &sm_mappings[store2]\n",
        "                  );\n",
        "\n",
        "                  if (!skip) calculateBin (\n",
        "                        bin_val2.y,\n",
        "                        &sm_mappings[store2 + blockDimx]\n",
        "                  );\n",
        "\n",
        "                  load_bin2 += input_pitch;\n",
        "                  store2 += width;\n",
        "                }\n",
        "        }\n",
        "\n",
        "        /*\n",
        "        if (does2) {\n",
        "          #pragma unroll\n",
        "          for (int i = 0; i < UNROLL; i++) {\n",
        "            uint2 bin_val2 = *load_bin2;\n",
        "\n",
        "            calculateBin (\n",
        "                bin_val2.x,\n",
        "                &sm_mappings[store2]\n",
        "            );\n",
        "\n",
        "            if (!skip) calculateBin (\n",
        "                           bin_val2.y,\n",
        "                           &sm_mappings[store2 + blockDimx]\n",
        "                        );\n",
        "\n",
        "            load_bin2 += input_pitch;\n",
        "            store2 += width;\n",
        "          }\n",
        "       }\n",
        "     */\n",
        "\n",
        "}\n",
        "\n"
};


//=========================================main=====================================//
#define DEFAULT_BLOCK_X         14
#define DEFAULT_PRESCAN_THREADS 512
#define DEFAULT_PRESCAN_BLOCKS_X    64
#define DEFAULT_FINAL_THREADS 512
#define UNROLL 16
#define UINT8_MAX 255
const int numIterations = 50;

//int runCPU(char *input) {
int runCPU(const char *input, double *exeTime) {
    unsigned int img_width, img_height;
    unsigned int histo_width, histo_height;

    flushed_printf("\tPreparing data: ");
    flushed_printf(input);
    FILE *f = fopen(input, "rb");
    int result = 0;
    result += fread(&img_width, sizeof(unsigned int), 1, f);
    result += fread(&img_height, sizeof(unsigned int), 1, f);
    result += fread(&histo_width, sizeof(unsigned int), 1, f);
    result += fread(&histo_height, sizeof(unsigned int), 1, f);
    if (result != 4) {
        fputs("Error reading input and output dimensions from file\n", stderr);
        return -1;
    }
    auto *img = (unsigned int *) malloc(img_width * img_height * sizeof(unsigned int));
    auto *histo = (unsigned char *) calloc(histo_width * histo_height, sizeof(unsigned char));
    result = fread(img, sizeof(unsigned int), img_width * img_height, f);
    fclose(f);

    if (result != img_width * img_height) {
        fputs("Error reading input array from file\n", stderr);
        return -1;
    }

    flushed_printf("\n\tImg size: %d x %d = %d", img_width, img_height, img_width * img_height);

    flushed_printf("\n\tHisto running...\n");
    double timer = getCurrentTime();
    int iter;
    for (iter = 0; iter < numIterations; iter++) {
        memset(histo, 0, histo_height * histo_width * sizeof(unsigned char));
        unsigned int i;

#pragma omp parallel for
        for (i = 0; i < img_width * img_height; ++i) {
            const unsigned int value = img[i];

#pragma omp critical
            if (histo[value] < UINT8_MAX) {
                ++histo[value];
            }
        }
    }

    timer = getCurrentTime() - timer;
    flushed_printf("\tHisto CPU done, time: %f sec\n", timer);
    *exeTime = timer;
    free(img);
    free(histo);
    return 0;
}

static const cl_uint PROGRAM_SOURCE_LEN = sizeof(PROGRAM_SOURCE) / sizeof(const char *);
cl_wrapper *cl;
cl_program program;
cl_kernel histo_prescan_kernel;
cl_kernel histo_intermediates_kernel;
cl_kernel histo_intermediates_kernel_compat;
cl_kernel histo_main_kernel;
cl_kernel histo_final_kernel;
cl_kernel chosenInterKernel;
cl_command_queue exeCmdQueue;

void initCL() {
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    exeCmdQueue = cl->createProfilingCmdQueue();
}

int runGPU(const char *input_str, double *exeTime) {
    initCL();
    unsigned int img_width, img_height;
    unsigned int histo_width, histo_height;
    unsigned int lmemKB;
    unsigned int prescanThreads = DEFAULT_PRESCAN_THREADS;
    unsigned int prescanBlockX = DEFAULT_PRESCAN_BLOCKS_X;
    unsigned int blockX = DEFAULT_BLOCK_X;
    unsigned int nThreads;
    unsigned int finalThreads = DEFAULT_FINAL_THREADS;
    unsigned int bins_per_block;

    FILE *f = fopen(input_str, "rb");
    int result = 0;

    result += fread(&img_width, sizeof(unsigned int), 1, f);
    result += fread(&img_height, sizeof(unsigned int), 1, f);
    result += fread(&histo_width, sizeof(unsigned int), 1, f);
    result += fread(&histo_height, sizeof(unsigned int), 1, f);

    if (result != 4) {
        fputs("Error reading input and output dimensions from file\n", stderr);
        return -1;
    }
//    unsigned int *img = (unsigned int *) malloc(img_width * img_height * sizeof(unsigned int));
//    unsigned char *histo = (unsigned char *) calloc(histo_width * histo_height, sizeof(unsigned char));
    ZeroCopyMem<unsigned int> img_z = init_zero_copy_region<unsigned int>(cl, img_width * img_height);
    ZeroCopyMem<unsigned char> histo_z = init_zero_copy_region<unsigned char>(cl, histo_width * histo_height);

//    result = fread(img, sizeof(unsigned int), img_width * img_height, f);
    result = fread(img_z.hostPtr, sizeof(unsigned int), img_width * img_height, f);
    fclose(f);
    if (result != img_width * img_height) {
        fputs("Error reading input array from file\n", stderr);
        return -1;
    }

    long unsigned int lmemSize = 0;
    clGetDeviceInfo(cl->device, CL_DEVICE_LOCAL_MEM_SIZE,
                    sizeof(cl_ulong), &lmemSize, nullptr);

    cl_uint workItemDimensions;
    clGetDeviceInfo(cl->device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
                    sizeof(cl_uint), &workItemDimensions, nullptr);

    size_t workItemSizes[workItemDimensions];
    clGetDeviceInfo(cl->device, CL_DEVICE_MAX_WORK_ITEM_SIZES,
                    workItemDimensions * sizeof(size_t), workItemSizes, nullptr);

    if (lmemSize >= 48 * 1024) lmemKB = 48;
    else if (lmemSize >= 24 * 1024) lmemKB = 24;
    else lmemKB = 8;

    bins_per_block = lmemKB * 1024;

    switch (lmemKB) {
        case 48:
            nThreads = 1024;
            break;
        case 24:
            nThreads = 768;
            break;
        default:
            nThreads = 512;
            break;
    }

    if ((workItemSizes[0] < 512) && (workItemSizes[0] >= 256)) {
        prescanThreads = DEFAULT_PRESCAN_THREADS / 2;
        prescanBlockX = DEFAULT_PRESCAN_BLOCKS_X * 2;
        blockX = DEFAULT_BLOCK_X * 2;
        nThreads = 256;
        finalThreads = DEFAULT_FINAL_THREADS / 2;
    }

    char compileOptions[1024];
    sprintf(compileOptions, " -D PRESCAN_THREADS=%u -D KB=%u -D UNROLL=%u -D BINS_PER_BLOCK=%u -D BLOCK_X=%u",
            prescanThreads, lmemKB, UNROLL, bins_per_block, blockX);
    program = cl->createProgramWithOptions(PROGRAM_SOURCE, PROGRAM_SOURCE_LEN, compileOptions);
    histo_prescan_kernel = cl->createKernel("histo_prescan_kernel", program);
    histo_intermediates_kernel = cl->createKernel(program, "histo_intermediates_kernel");
    histo_intermediates_kernel_compat = cl->createKernel(program, "histo_intermediates_kernel_compat");
    histo_main_kernel = cl->createKernel(program, "histo_main_kernel");
    histo_final_kernel = cl->createKernel(program, "histo_final_kernel");


    int even_width = ((img_width + 1) / 2) * 2;

    ZeroCopyMem<unsigned int> input_z = init_zero_copy_region<unsigned int>
            (cl, even_width * (((img_height + UNROLL - 1) / UNROLL) * UNROLL));
    ZeroCopyMem<unsigned int> ranges_z = init_zero_copy_region<unsigned int>(cl, 2);
    ZeroCopyMem<unsigned char> sm_mappings_z = init_zero_copy_region<unsigned char>(cl, img_width * img_height * 4);
    ZeroCopyMem<unsigned int> global_subhisto_z = init_zero_copy_region<unsigned int>(cl, img_width * histo_height);
    ZeroCopyMem<unsigned short> global_histo_z = init_zero_copy_region<unsigned short>(cl, img_width * histo_height);
    ZeroCopyMem<unsigned int> global_overflow_z = init_zero_copy_region<unsigned int>(cl, img_width * histo_height);
    ZeroCopyMem<unsigned char> final_histo_z = init_zero_copy_region<unsigned char>(cl, img_width * histo_height);

    unsigned int *zeroData;
    zeroData = (unsigned int *) calloc(img_width * histo_height, sizeof(unsigned int));
    if (zeroData == NULL) {
        fprintf(stderr, "Failed to allocate %ld bytes of memory on host!\n",
                sizeof(unsigned int) * img_width * histo_height);
        exit(1);
    }

    for (int y = 0; y < img_height; y++) {
        memcpy(&input_z.hostPtr[y * even_width], &img_z.hostPtr[y * img_width], img_width * sizeof(unsigned int));
//        cl_wrapper::checkError(clEnqueueWriteBuffer(clCommandQueue, input, CL_TRUE,
//                                              y * even_width * sizeof(unsigned int), // Offset in bytes
//                                              img_width * sizeof(unsigned int), // Size of data to write
//                                              &img[y * img_width], // Host Source
//                                              0, NULL, NULL));
    }

    double time = getCurrentTime();
    unsigned int img_dim = img_height * img_width;
    clSetKernelArg(histo_prescan_kernel, 0, sizeof(cl_mem), &input_z.deviceBuffer);
    clSetKernelArg(histo_prescan_kernel, 1, sizeof(unsigned int), &img_dim);
    clSetKernelArg(histo_prescan_kernel, 2, sizeof(cl_mem), (void *) &ranges_z.deviceBuffer);

    unsigned int half_width = (img_width + 1) / 2;
    clSetKernelArg(histo_intermediates_kernel, 0, sizeof(cl_mem), (void *) &input_z.deviceBuffer);
    clSetKernelArg(histo_intermediates_kernel, 1, sizeof(unsigned int), &img_height);
    clSetKernelArg(histo_intermediates_kernel, 2, sizeof(unsigned int), &img_width);
    clSetKernelArg(histo_intermediates_kernel, 3, sizeof(unsigned int), &half_width);
    clSetKernelArg(histo_intermediates_kernel, 4, sizeof(cl_mem), (void *) &sm_mappings_z.deviceBuffer);

    clSetKernelArg(histo_intermediates_kernel_compat, 0, sizeof(cl_mem), (void *) &input_z.deviceBuffer);
    clSetKernelArg(histo_intermediates_kernel_compat, 1, sizeof(unsigned int), &img_height);
    clSetKernelArg(histo_intermediates_kernel_compat, 2, sizeof(unsigned int), &img_width);
    clSetKernelArg(histo_intermediates_kernel_compat, 3, sizeof(unsigned int), &half_width);
    clSetKernelArg(histo_intermediates_kernel_compat, 4, sizeof(cl_mem), (void *) &sm_mappings_z.deviceBuffer);


    clSetKernelArg(histo_main_kernel, 0, sizeof(cl_mem), (void *) &sm_mappings_z.deviceBuffer);
    clSetKernelArg(histo_main_kernel, 1, sizeof(unsigned int), &img_dim);

    clSetKernelArg(histo_main_kernel, 4, sizeof(unsigned int), &histo_height);
    clSetKernelArg(histo_main_kernel, 5, sizeof(unsigned int), &histo_width);
    clSetKernelArg(histo_main_kernel, 6, sizeof(cl_mem), (void *) &global_subhisto_z.deviceBuffer);
    clSetKernelArg(histo_main_kernel, 7, sizeof(cl_mem), (void *) &global_histo_z.deviceBuffer);
    clSetKernelArg(histo_main_kernel, 8, sizeof(cl_mem), (void *) &global_overflow_z.deviceBuffer);


    clSetKernelArg(histo_final_kernel, 2, sizeof(unsigned int), &histo_height);
    clSetKernelArg(histo_final_kernel, 3, sizeof(unsigned int), &histo_width);
    clSetKernelArg(histo_final_kernel, 4, sizeof(cl_mem), (void *) &global_subhisto_z.deviceBuffer);
    clSetKernelArg(histo_final_kernel, 5, sizeof(cl_mem), (void *) &global_histo_z.deviceBuffer);
    clSetKernelArg(histo_final_kernel, 6, sizeof(cl_mem), (void *) &global_overflow_z.deviceBuffer);
    clSetKernelArg(histo_final_kernel, 7, sizeof(cl_mem), (void *) &final_histo_z.deviceBuffer);

    size_t prescan_localWS[1] = {prescanThreads};
    size_t prescan_globalWS[1] = {prescanBlockX * prescan_localWS[0]};

    size_t inter_localWS[1] = {half_width};
    size_t inter_globalWS[1] = {((img_height + UNROLL - 1) / UNROLL) * inter_localWS[0]};

    size_t main_localWS[2] = {nThreads, 1};
    size_t main_globalWS[2] = {blockX * main_localWS[0], 0};

    size_t final_localWS[1] = {finalThreads};
    size_t final_globalWS[1] = {blockX * 3 * final_localWS[0]};

    if (half_width > workItemSizes[0]) {
        chosenInterKernel = histo_intermediates_kernel_compat;
        inter_localWS[0] = workItemSizes[0];
        inter_globalWS[0] = ((img_height + UNROLL - 1) / UNROLL) * inter_localWS[0];
    } else {
        chosenInterKernel = histo_intermediates_kernel;
    }

//    cl_event event[4];
    for (int iter = 0; iter < numIterations; iter++) {
        unsigned int ranges_h[2] = {UINT32_MAX, 0};
        memcpy(&ranges_z.hostPtr[0], &ranges_h[0], 2 * sizeof(unsigned int));
        cl_wrapper::checkError(clEnqueueNDRangeKernel(
//                exeCmdQueue, histo_prescan_kernel, 1, 0, prescan_globalWS, prescan_localWS, 0, 0, &event[0]));
                exeCmdQueue, histo_prescan_kernel, 1, 0, prescan_globalWS, prescan_localWS, 0, 0, 0));
//        clWaitForEvents(1, &event[0]);
        memcpy(&ranges_h[0], &ranges_z.hostPtr[0], 2 * sizeof(unsigned int));
        memcpy(&global_subhisto_z.hostPtr[0], zeroData, img_width * histo_height * sizeof(unsigned int));
        cl_wrapper::checkError(
                clEnqueueNDRangeKernel(exeCmdQueue, chosenInterKernel, 1, 0, inter_globalWS, inter_localWS, 0, 0,
//                                       &event[1]));
                                       0));
//        clWaitForEvents(1, &event[1]);

        main_globalWS[1] = ranges_h[1] - ranges_h[0] + 1;
        cl_wrapper::checkError(clSetKernelArg(histo_main_kernel, 2, sizeof(unsigned int), &ranges_h[0]));
        cl_wrapper::checkError(clSetKernelArg(histo_main_kernel, 3, sizeof(unsigned int), &ranges_h[1]));
        cl_wrapper::checkError(clSetKernelArg(histo_final_kernel, 0, sizeof(unsigned int), &ranges_h[0]));
        cl_wrapper::checkError(clSetKernelArg(histo_final_kernel, 1, sizeof(unsigned int), &ranges_h[1]));

        cl_wrapper::checkError(
                clEnqueueNDRangeKernel(exeCmdQueue, histo_main_kernel, 2, 0, main_globalWS, main_localWS, 0, 0,
//                                       &event[2]));
                                       0));
//        clWaitForEvents(1, &event[2]);

        cl_wrapper::checkError(clEnqueueNDRangeKernel(
//                exeCmdQueue, histo_final_kernel, 1, 0, final_globalWS, final_localWS, 0, 0, &event[3]));
                exeCmdQueue, histo_final_kernel, 1, 0, final_globalWS, final_localWS, 0, 0, 0));
//        clWaitForEvents(1, &event[3]);
//        printf("%d\n", iter);
    }
//    for (auto &i : event) {
//        clReleaseEvent(i);
//    }

    time = getCurrentTime() - time;
    flushed_printf("\thisto GPU done, time: %f sec\n", time);

    *exeTime = time;

    cl_wrapper::checkError(clReleaseMemObject(input_z.deviceBuffer));
    cl_wrapper::checkError(clReleaseMemObject(ranges_z.deviceBuffer));
    cl_wrapper::checkError(clReleaseMemObject(sm_mappings_z.deviceBuffer));
    cl_wrapper::checkError(clReleaseMemObject(global_subhisto_z.deviceBuffer));
    cl_wrapper::checkError(clReleaseMemObject(global_histo_z.deviceBuffer));
    cl_wrapper::checkError(clReleaseMemObject(global_overflow_z.deviceBuffer));
    cl_wrapper::checkError(clReleaseMemObject(final_histo_z.deviceBuffer));
    cl_wrapper::checkError(clReleaseMemObject(img_z.deviceBuffer));
    cl_wrapper::checkError(clReleaseMemObject(histo_z.deviceBuffer));
    free(zeroData);
    return 0;
}

int main(int argc, char **argv) {
    char const *input = "datasets/histo/large/input/img.bin";
    char const *input_small = "datasets/histo/default/input/img.bin";
//    double test;
//    runGPU(input, &test);
    benchmark(input, input_small, runCPU, runGPU);
    return 0;
}




















