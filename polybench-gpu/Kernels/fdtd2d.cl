MSTRINGIFY(
\n/**
\n * fdtd2d.cl: This file is part of the PolyBench/GPU 1.0 test suite.
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
\n__kernel void fdtd_kernel1(__global DATA_TYPE *_fict_, __global DATA_TYPE *ex, __global DATA_TYPE *ey, __global DATA_TYPE *hz, int t, int nx, int ny) 
\n{    
\n    int j = get_global_id(0);
\n    int i = get_global_id(1);
\n
\n    if ((i < nx) && (j < ny))
\n    {
\n        int tid = i * ny + j;
\n
\n        if (i == 0) 
\n        {
\n            ey[i * ny + j] = _fict_[t];
\n        }
\n        else
\n        { 
\n            ey[i * ny + j] = ey[i * ny + j] - 0.5*(hz[i * ny + j] - hz[(i-1) * ny + j]);
\n        }
\n    }
\n}
\n
\n
\n__kernel void fdtd_kernel2(__global DATA_TYPE *ex, __global DATA_TYPE *ey, __global DATA_TYPE *hz, int nx, int ny) 
\n{    
\n    int j = get_global_id(0);
\n    int i = get_global_id(1);
\n
\n    if ((i < nx) && (j < ny) && (j > 0))
\n    {
\n        ex[i * (ny+1) + j] = ex[i * (ny+1) + j] - 0.5*(hz[i * ny + j] - hz[i * ny + (j-1)]);
\n    }
\n}
\n
\n
\n__kernel void fdtd_kernel3(__global DATA_TYPE *ex, __global DATA_TYPE *ey, __global DATA_TYPE *hz, int nx, int ny) 
\n{    
\n    int j = get_global_id(0);
\n    int i = get_global_id(1);
\n
\n    if ((i < nx) && (j < ny))
\n    {
\n        hz[i * ny + j] = hz[i * ny + j] - 0.7*(ex[i * (ny+1) + (j+1)] - ex[i * (ny+1) + j] + ey[(i + 1) * ny + j] - ey[i * ny + j]);
\n    }
\n} 
);
