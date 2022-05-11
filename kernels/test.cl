
__kernel void kernel1D(__global float* record){
    int i = get_global_id(0);
    record[i] = (float) i/100;
}

__kernel void kernel2D(__global float* record){
    int i = get_global_id(0);
    int j = get_global_id(1);
    int size_i = get_global_size(0);
    int loc = i + j*size_i;
    record[loc] = (float) loc/100;
}

__kernel void kernel3D(__global float* record){
    int i = get_global_id(0);
    int j = get_global_id(1);
    int k = get_global_id(2);
    int size_i = get_global_size(0);
    int size_j = get_global_size(1);
    int loc = i + j*size_i + k*size_j;
    record[loc] = loc/1000;
}
