//
// Created by pfxu on 20/8/2019.
//

#ifndef MOBILEHETEROGENOUSPROJECT_CLION_MRI_GRIDDING_H
#define MOBILEHETEROGENOUSPROJECT_CLION_MRI_GRIDDING_H

#define PI 3.14159265
typedef struct {
    int numSamples;
    int aquisitionMatrixSize[3];
    int reconstructionMatrixSize[3];
    float kMax[3];
    int gridSize[3];
    float oversample;
    float kernelWidth;
    int binsize;
    int useLUT;
} parameters;

typedef struct {
    float real;
    float imag;
    float kX;
    float kY;
    float kZ;
    float sdc;
} ReconstructionSample;

typedef struct {
    float real;
    float imag;
} cmplx;

#define max(x, y) ((x<y)?y:x)
#define min(x, y) ((x>y)?y:x)

#define PI 3.14159265359

float kernel_value_CPU(float v) {

    float rValue = 0;

    const float z = v * v;

    // polynomials taken from http://ccrma.stanford.edu/CCRMA/Courses/422/projects/kbd/kbdwindow.cpp
    float num = (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * 0.210580722890567e-22f +
                                                                                  0.380715242345326e-19f) +
                                                                             0.479440257548300e-16f) +
                                                                        0.435125971262668e-13f) +
                                                                   0.300931127112960e-10f) + 0.160224679395361e-7f) +
                                                         0.654858370096785e-5f) + 0.202591084143397e-2f) +
                                               0.463076284721000e0f) + 0.754337328948189e2f) + 0.830792541809429e4f) +
                                0.571661130563785e6f) + 0.216415572361227e8f) + 0.356644482244025e9f) +
                 0.144048298227235e10f);
    float den = (z * (z * (z - 0.307646912682801e4f) + 0.347626332405882e7f) - 0.144048298227235e10f);
    rValue = -num / den;
    return rValue;
}

void calculateLUT(float beta, float width, float **LUT, unsigned int *sizeLUT) {
    float v;
    float cutoff2 = (width * width) / 4.0;

    unsigned int size;

    if (width > 0) {
        // compute size of LUT based on kernel width
        size = (unsigned int) (10000 * width);

        // allocate memory
        (*LUT) = (float *) malloc(size * sizeof(float));

        unsigned int k;

#pragma omp parallel for private(v)
        for (k = 0; k < size; ++k) {
            // compute value to evaluate kernel at
            // v in the range 0:(_width/2)^2
            v = (((float) k) / ((float) size)) * cutoff2;

            // compute kernel value and store
            (*LUT)[k] = kernel_value_CPU(beta * sqrt(1.0 - (v / cutoff2)));
        }
        (*sizeLUT) = size;
    }
}

float kernel_value_LUT(float v, float *LUT, int sizeLUT, float _1overCutoff2) {
    unsigned int k0;
    float v0;

    v *= (float) sizeLUT;
    k0 = (unsigned int) (v * _1overCutoff2);
    v0 = ((float) k0) / _1overCutoff2;
    return LUT[k0] + ((v - v0) * (LUT[k0 + 1] - LUT[k0]) / _1overCutoff2);
}

int gridding_Gold(unsigned int n, parameters params, ReconstructionSample *sample, float *LUT, unsigned int sizeLUT,
                  cmplx *gridData, float *sampleDensity) {

    unsigned int NxL, NxH;
    unsigned int NyL, NyH;
    unsigned int NzL, NzH;

    int nx;
    int ny;
    int nz;

    float w;
    unsigned int idx;
    unsigned int idx0;

    unsigned int idxZ;
    unsigned int idxY;

    float Dx2[100];
    float Dy2[100];
    float Dz2[100];
    float *dx2 = nullptr;
    float *dy2 = nullptr;
    float *dz2 = nullptr;

    float dy2dz2;
    float v;

    unsigned int size_x = params.gridSize[0];
    unsigned int size_y = params.gridSize[1];
    unsigned int size_z = params.gridSize[2];

    float cutoff = ((float) (params.kernelWidth)) / 2.0; // cutoff radius
    float cutoff2 = cutoff * cutoff;                    // square of cutoff radius
    float _1overCutoff2 = 1 / cutoff2;                  // 1 over square of cutoff radius

    float beta = PI * sqrt(4 * params.kernelWidth * params.kernelWidth / (params.oversample * params.oversample) *
                           (params.oversample - .5) * (params.oversample - .5) - .8);

    int i;

#pragma omp parallel for private(NxL, NxH, NyL, NyH, NzL, NzH, dz2, nz, dx2, \
                 nx, dy2, ny, idxZ, idxY, dy2dz2, idx0, v, idx, w)

    for (i = 0; i < n; i++) {
        ReconstructionSample pt = sample[i];

        float kx = pt.kX;
        float ky = pt.kY;
        float kz = pt.kZ;

        NxL = max((kx - cutoff), 0.0);
        NxH = min((kx + cutoff), size_x - 1.0);

        NyL = max((ky - cutoff), 0.0);
        NyH = min((ky + cutoff), size_y - 1.0);

        NzL = max((kz - cutoff), 0.0);
        NzH = min((kz + cutoff), size_z - 1.0);

        if ((pt.real != 0.0 || pt.imag != 0.0) && pt.sdc != 0.0) {
            for (dz2 = Dz2, nz = NzL; nz <= NzH; ++nz, ++dz2) {
                *dz2 = ((kz - nz) * (kz - nz));
            }
            for (dx2 = Dx2, nx = NxL; nx <= NxH; ++nx, ++dx2) {
                *dx2 = ((kx - nx) * (kx - nx));
            }
            for (dy2 = Dy2, ny = NyL; ny <= NyH; ++ny, ++dy2) {
                *dy2 = ((ky - ny) * (ky - ny));
            }

            idxZ = (NzL - 1) * size_x * size_y;
            for (dz2 = Dz2, nz = NzL; nz <= NzH; ++nz, ++dz2) {
                /* linear offset into 3-D matrix to get to zposition */
                idxZ += size_x * size_y;

                idxY = (NyL - 1) * size_x;

                /* loop over x indexes, but only if curent distance is close enough (distance will increase by adding x&y distance) */
                if ((*dz2) < cutoff2) {
                    for (dy2 = Dy2, ny = NyL; ny <= NyH; ++ny, ++dy2) {
                        /* linear offset IN ADDITION to idxZ to get to Y position */
                        idxY += size_x;
                        dy2dz2 = (*dz2) + (*dy2);
                        idx0 = idxY + idxZ;
                        /* loop over y indexes, but only if curent distance is close enough (distance will increase by adding y distance) */
                        if (dy2dz2 < cutoff2) {
                            for (dx2 = Dx2, nx = NxL; nx <= NxH; ++nx, ++dx2) {
                                /* value to evaluate kernel at */
                                v = dy2dz2 + (*dx2);
                                if (v < cutoff2) {
                                    /* linear index of (x,y,z) point */
                                    idx = nx + idx0;
                                    /* kernel weighting value */
                                    if (params.useLUT) {
                                        w = kernel_value_LUT(v, LUT, sizeLUT, _1overCutoff2) * pt.sdc;
                                    } else {
                                        w = kernel_value_CPU(beta * sqrt(1.0 - (v * _1overCutoff2))) * pt.sdc;
                                    }
                                    /* grid data */
#pragma omp critical (c1)
                                    gridData[idx].real += (w * pt.real);
#pragma omp critical (c2)
                                    gridData[idx].imag += (w * pt.imag);
                                    /* estimate sample density */
#pragma omp critical (c3)
                                    sampleDensity[idx] += 1.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

void setParameters(FILE *file, parameters *p) {
    fscanf(file, "aquisition.numsamples=%d\n", &(p->numSamples));
    fscanf(file, "aquisition.kmax=%f %f %f\n", &(p->kMax[0]), &(p->kMax[1]), &(p->kMax[2]));
    fscanf(file, "aquisition.matrixSize=%d %d %d\n", &(p->aquisitionMatrixSize[0]), &(p->aquisitionMatrixSize[1]),
           &(p->aquisitionMatrixSize[2]));
    fscanf(file, "reconstruction.matrixSize=%d %d %d\n", &(p->reconstructionMatrixSize[0]),
           &(p->reconstructionMatrixSize[1]), &(p->reconstructionMatrixSize[2]));
    fscanf(file, "gridding.matrixSize=%d %d %d\n", &(p->gridSize[0]), &(p->gridSize[1]), &(p->gridSize[2]));
    fscanf(file, "gridding.oversampling=%f\n", &(p->oversample));
    fscanf(file, "kernel.width=%f\n", &(p->kernelWidth));
    fscanf(file, "kernel.useLUT=%d\n", &(p->useLUT));

    printf("  Number of samples = %d\n", p->numSamples);
    printf("  Grid Size = %dx%dx%d\n", p->gridSize[0], p->gridSize[1], p->gridSize[2]);
    printf("  Input Matrix Size = %dx%dx%d\n", p->aquisitionMatrixSize[0], p->aquisitionMatrixSize[1],
           p->aquisitionMatrixSize[2]);
    printf("  Recon Matrix Size = %dx%dx%d\n", p->reconstructionMatrixSize[0], p->reconstructionMatrixSize[1],
           p->reconstructionMatrixSize[2]);
    printf("  Kernel Width = %f\n", p->kernelWidth);
    printf("  KMax = %.2f %.2f %.2f\n", p->kMax[0], p->kMax[1], p->kMax[2]);
    printf("  Oversampling = %f\n", p->oversample);
    printf("  GPU Binsize = %d\n", p->binsize);
    printf("  Use LUT = %s\n", (p->useLUT) ? "Yes" : "No");
}


unsigned int readSampleData(parameters params, FILE *uksdata_f, ReconstructionSample *samples) {
    unsigned int i;

    for (i = 0; i < params.numSamples; i++) {
        if (feof(uksdata_f)) {
            break;
        }
        fread((void *) &(samples[i]), sizeof(ReconstructionSample), 1, uksdata_f);
    }

    float kScale[3];
    kScale[0] = (float) (params.aquisitionMatrixSize[0]) /
                ((float) (params.reconstructionMatrixSize[0]) * (float) (params.kMax[0]));
    kScale[1] = (float) (params.aquisitionMatrixSize[1]) /
                ((float) (params.reconstructionMatrixSize[1]) * (float) (params.kMax[1]));
    kScale[2] = (float) (params.aquisitionMatrixSize[2]) /
                ((float) (params.reconstructionMatrixSize[2]) * (float) (params.kMax[2]));

    int size_x = params.gridSize[0];
    int size_y = params.gridSize[1];
    int size_z = params.gridSize[2];

    float ax = (kScale[0] * (size_x - 1)) / 2.0;
    float bx = (float) (size_x - 1) / 2.0;

    float ay = (kScale[1] * (size_y - 1)) / 2.0;
    float by = (float) (size_y - 1) / 2.0;

    float az = (kScale[2] * (size_z - 1)) / 2.0;
    float bz = (float) (size_z - 1) / 2.0;

    int n;
    for (n = 0; n < i; n++) {
        samples[n].kX = floor((samples[n].kX * ax) + bx);
        samples[n].kY = floor((samples[n].kY * ay) + by);
        samples[n].kZ = floor((samples[n].kZ * az) + bz);
    }

    return i;
}

static const char *PROGRAM_SOURCE[] = {
        "#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable\n",
        "\n",
        "typedef struct{\n",
        "  float real;\n",
        "  float imag;\n",
        "  float kX;\n",
        "  float kY;\n",
        "  float kZ;\n",
        "  float sdc;\n",
        "} ReconstructionSample;\n",
        "\n",
        "#define TILE 64\n",
        "#define LOG_TILE 6\n",
        "\n",
        "__kernel void binning_kernel (unsigned int n,\n",
        "  __global ReconstructionSample* sample_g,\n",
        "  __global unsigned int* idxKey_g,\n",
        "  __global unsigned int* idxValue_g,\n",
        "  __global unsigned int* binCount_g,\n",
        "  unsigned int binsize, unsigned int gridNumElems){\n",
        "  unsigned int key;\n",
        "  unsigned int sampleIdx = get_global_id(0);\n",
        "  ReconstructionSample pt;\n",
        "  unsigned int binIdx;\n",
        "  unsigned int count;\n",
        "  if (sampleIdx < n){\n",
        "    pt = sample_g[sampleIdx];\n",
        "    binIdx = (unsigned int)(pt.kZ)*((int) ( SIZE_XY_VAL )) + (unsigned int)(pt.kY)*((int)( GRIDSIZE_VAL1 )) + (unsigned int)(pt.kX);\n",
        "    count = atom_add(binCount_g+binIdx, 1);\n",
        "    if (count < binsize){\n",
        "      key = binIdx;\n",
        "    } else {\n",
        "      atom_sub(binCount_g+binIdx, 1);\n",
        "      key = gridNumElems;\n",
        "    }\n",
        "    idxKey_g[sampleIdx] = key;\n",
        "    idxValue_g[sampleIdx] = sampleIdx;\n",
        "  }\n",
        "}\n",
        "\n",
        "__kernel void reorder_kernel(int n,\n",
        " __global unsigned int* idxValue_g,\n",
        " __global ReconstructionSample* samples_g,\n",
        " __global ReconstructionSample* sortedSample_g){\n",
        "  unsigned int index = get_global_id(0);\n",
        "  unsigned int old_index;\n",
        "  ReconstructionSample pt;\n",
        "  if (index < n){\n",
        "    old_index = idxValue_g[index];\n",
        "    pt = samples_g[old_index];\n",
        "    sortedSample_g[index] = pt;\n",
        "  }\n",
        "}\n",
        "float kernel_value(float v){\n",
        "  float rValue = 0;\n",
        "  float z = v*v;\n",
        "  float num = (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z*\n",
        "    (z* 0.210580722890567e-22f  + 0.380715242345326e-19f ) +\n",
        "    0.479440257548300e-16f) + 0.435125971262668e-13f ) +\n",
        "  0.300931127112960e-10f) + 0.160224679395361e-7f  ) +\n",
        "  0.654858370096785e-5f)  + 0.202591084143397e-2f  ) +\n",
        "  0.463076284721000e0f)   + 0.754337328948189e2f   ) +\n",
        "  0.830792541809429e4f)   + 0.571661130563785e6f   ) +\n",
        "  0.216415572361227e8f)   + 0.356644482244025e9f   ) +\n",
        "  0.144048298227235e10f);\n",
        "\n",
        "  float den = (z*(z*(z-0.307646912682801e4f)+0.347626332405882e7f)-0.144048298227235e10f);\n",
        "  rValue = native_divide(-num,den);\n",
        "  return rValue;\n",
        "}\n",
        "\n",
        "__kernel void gridding_GPU (__global ReconstructionSample* sample_g,\n",
        "  __global unsigned int* binStartAddr_g,\n",
        "  __global float2* gridData_g,\n",
        "  __global float* sampleDensity_g,\n",
        "  float beta\n",
        "  ){\n",
        "  __local ReconstructionSample sharedBin[TILE];\n",
        "  const int flatIdx = get_local_id(2)*get_local_size(1)*get_local_size(0) + get_local_id(1)*get_local_size(0) + get_local_id(0);\n",
        "  const int z0 = get_local_size(2)*(get_group_id(1)/(GRIDSIZE_VAL2/get_local_size(1)));\n",
        "  const int y0 = get_local_size(1)*(get_group_id(1)%(GRIDSIZE_VAL2/get_local_size(1)));\n",
        "  const int x0 = get_group_id(0)*get_local_size(0);\n",
        "  const int X  = x0+get_local_id(0);\n",
        "  const int Y  = y0+get_local_id(1);\n",
        "  const int Z  = z0+get_local_id(2);\n",
        "  const int xl = x0-CEIL_CUTOFF_VAL;\n",
        "  const int xL = (xl < 0) ? 0 : xl;\n",
        "  const int xh = x0+get_local_size(0)+CUTOFF_VAL;\n",
        "  const int xH = (xh >= GRIDSIZE_VAL1) ? GRIDSIZE_VAL1-1 : xh;\n",
        "  const int yl = y0-CEIL_CUTOFF_VAL;\n",
        "  const int yL = (yl < 0) ? 0 : yl;\n",
        "  const int yh = y0+get_local_size(1)+CUTOFF_VAL;\n",
        "  const int yH = (yh >= GRIDSIZE_VAL2) ? GRIDSIZE_VAL2-1 : yh;\n",
        "  const int zl = z0-CEIL_CUTOFF_VAL;\n",
        "  const int zL = (zl < 0) ? 0 : zl;\n",
        "  const int zh = z0+get_local_size(2)+CUTOFF_VAL;\n",
        "  const int zH = (zh >= GRIDSIZE_VAL3) ? GRIDSIZE_VAL3-1 : zh;\n",
        "  const int idx = Z*SIZE_XY_VAL + Y*GRIDSIZE_VAL1 + X;\n",
        "  float2 pt = (float2) (0.0f, 0.0f);\n",
        "  float density = 0.0f;\n",
        "  for (int z = zL; z <= zH; z++){\n",
        "    for (int y = yL; y <= yH; y++){\n",
        "      __global const unsigned int *addr = binStartAddr_g+z*SIZE_XY_VAL+ y*GRIDSIZE_VAL1;\n",
        "      const unsigned int start = *(addr+xL);\n",
        "      const unsigned int end   = *(addr+xH+1);\n",
        "      const unsigned int delta = end-start;\n",
        "      for (int x = 0; x < ((delta+TILE-1)>>LOG_TILE); x++){\n",
        "        int tileSize = ((delta-(x<<LOG_TILE)) > TILE) ? TILE : (delta-(x<<LOG_TILE));\n",
        "        int globalIdx = flatIdx+(x<<LOG_TILE);\n",
        "        barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "        if(flatIdx < tileSize){\n",
        "          sharedBin[flatIdx] = sample_g[start+globalIdx];\n",
        "        }\n",
        "        barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "        for (int j=0; j< tileSize; j++){\n",
        "          const float real = sharedBin[j].real;\n",
        "          const float imag = sharedBin[j].imag;\n",
        "          const float sdc = sharedBin[j].sdc;\n",
        "          if((real != 0.0f || imag != 0.0f) && sdc != 0.0f){\n",
        "            float v = (sharedBin[j].kX-X)*(sharedBin[j].kX-X);\n",
        "            v += (sharedBin[j].kY-Y)*(sharedBin[j].kY-Y);\n",
        "            v += (sharedBin[j].kZ-Z)*(sharedBin[j].kZ-Z);\n",
        "            if(v<CUTOFF2_VAL){\n",
        "             const float w = kernel_value(beta*sqrt(1.0f-(v*ONE_OVER_CUTOFF2_VAL))) *sdc;\n",
        "             pt.x += w*real;\n",
        "             pt.y += w*imag;\n",
        "             density += 1.0f;\n",
        "           }\n",
        "         }\n",
        "       }\n",
        "     }\n",
        "   }\n",
        " }\n",
        " gridData_g[idx] = pt;\n",
        " sampleDensity_g[idx] = density;\n",
        "}\n",
        "\n",
        "#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable\n",
        "\n",
        "#define UINT32_MAX 4294967295\n",
        "#define BITS 4\n",
        "#define LNB 4\n",
        "\n",
        "#define SORT_BS 256\n",
        "\n",
        "//#define CONFLICT_FREE_OFFSET(index) ((index) >> LNB + (index) >> (2*LNB))\n",
        "#define CONFLICT_FREE_OFFSET(index) (((unsigned int)(index) >> min((unsigned int)((LNB)+(index)), (unsigned int)((32-(2*LNB)))))>>(2*LNB))\n",
        "#define BLOCK_P_OFFSET (4*SORT_BS+1+(4*SORT_BS+1)/16+(4*SORT_BS+1)/64)\n",
        "\n",
        "void scan (__local unsigned int s_data[BLOCK_P_OFFSET]){\n",
        "  unsigned int thid = get_local_id(0);\n",
        "  barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "  s_data[2*thid+1+CONFLICT_FREE_OFFSET(2*thid+1)] += s_data[2*thid+CONFLICT_FREE_OFFSET(2*thid)];\n",
        "  s_data[2*(get_local_size(0)+thid)+1+CONFLICT_FREE_OFFSET(2*(get_local_size(0)+thid)+1)] += s_data[2*(get_local_size(0)+thid)+CONFLICT_FREE_OFFSET(2*(get_local_size(0)+thid))];\n",
        "  unsigned int stride = 2;\n",
        "  for (unsigned int d = get_local_size(0); d > 0; d >>= 1)\n",
        "  {\n",
        "    barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "    if (thid < d)\n",
        "    {\n",
        "      unsigned int i  = 2*stride*thid;\n",
        "      unsigned int ai = i + stride - 1;\n",
        "      unsigned int bi = ai + stride;\n",
        "\n",
        "      ai += CONFLICT_FREE_OFFSET(ai);\n",
        "      bi += CONFLICT_FREE_OFFSET(bi);\n",
        "\n",
        "      s_data[bi] += s_data[ai];\n",
        "    }\n",
        "    stride *= 2;\n",
        "  }\n",
        "  if (thid == 0){\n",
        "    unsigned int last = 4*get_local_size(0)-1;\n",
        "    last += CONFLICT_FREE_OFFSET(last);\n",
        "    s_data[4*get_local_size(0)+CONFLICT_FREE_OFFSET(4*get_local_size(0))] = s_data[last];\n",
        "    s_data[last] = 0;\n",
        "  }\n",
        "\n",
        "  for (unsigned int d = 1; d <= get_local_size(0); d *= 2)\n",
        "  {\n",
        "    stride >>= 1;\n",
        "    barrier(CLK_LOCAL_MEM_FENCE );\n",
        "    if (thid < d)\n",
        "    {\n",
        "      unsigned int i  = 2*stride*thid;\n",
        "      unsigned int ai = i + stride - 1;\n",
        "      unsigned int bi = ai + stride;\n",
        "      ai += CONFLICT_FREE_OFFSET(ai);\n",
        "      bi += CONFLICT_FREE_OFFSET(bi);\n",
        "      unsigned int t  = s_data[ai];\n",
        "      s_data[ai] = s_data[bi];\n",
        "      s_data[bi] += t;\n",
        "    }\n",
        "  }\n",
        "  barrier(CLK_LOCAL_MEM_FENCE );\n",
        "  unsigned int temp = s_data[2*thid+CONFLICT_FREE_OFFSET(2*thid)];\n",
        "  s_data[2*thid+CONFLICT_FREE_OFFSET(2*thid)] = s_data[2*thid+1+CONFLICT_FREE_OFFSET(2*thid+1)];\n",
        "  s_data[2*thid+1+CONFLICT_FREE_OFFSET(2*thid+1)] += temp;\n",
        "  unsigned int temp2 = s_data[2*(get_local_size(0)+thid)+CONFLICT_FREE_OFFSET(2*(get_local_size(0)+thid))];\n",
        "  s_data[2*(get_local_size(0)+thid)+CONFLICT_FREE_OFFSET(2*(get_local_size(0)+thid))] = s_data[2*(get_local_size(0)+thid)+1+CONFLICT_FREE_OFFSET(2*(get_local_size(0)+thid)+1)];\n",
        "  s_data[2*(get_local_size(0)+thid)+1+CONFLICT_FREE_OFFSET(2*(get_local_size(0)+thid)+1)] += temp2;\n",
        "  barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "}\n",
        "\n",
        "__kernel void splitSort(int numElems, int iter,\n",
        " __global unsigned int* keys,\n",
        " __global unsigned int* values,\n",
        " __global unsigned int* histo){\n",
        "  __local unsigned int flags[BLOCK_P_OFFSET];\n",
        "  __local unsigned int histo_s[1<<BITS];\n",
        "  const unsigned int tid = get_local_id(0);\n",
        "  const unsigned int gid = get_group_id(0)*4*SORT_BS+4*get_local_id(0);\n",
        "  uint4 lkey = { UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX};\n",
        "  uint4 lvalue;\n",
        "  if (gid < numElems){\n",
        "    lkey = *((__global uint4*)(keys+gid));\n",
        "    lvalue = *((__global uint4*)(values+gid));\n",
        "  }\n",
        "  if(tid < (1<<BITS)){\n",
        "    histo_s[tid] = 0;\n",
        "  }\n",
        "  barrier(CLK_LOCAL_MEM_FENCE );\n",
        "  atom_add(histo_s+((lkey.x&((1<<(BITS*(iter+1)))-1))>>(BITS*iter)),1);\n",
        "  atom_add(histo_s+((lkey.y&((1<<(BITS*(iter+1)))-1))>>(BITS*iter)),1);\n",
        "  atom_add(histo_s+((lkey.z&((1<<(BITS*(iter+1)))-1))>>(BITS*iter)),1);\n",
        "  atom_add(histo_s+((lkey.w&((1<<(BITS*(iter+1)))-1))>>(BITS*iter)),1);\n",
        "  uint4 index = (uint4) (4*tid, 4*tid+1, 4*tid+2, 4*tid+3);\n",
        "  for (int i=BITS*iter; i<BITS*(iter+1);i++){\n",
        "    const uint4 flag = (uint4) ( (lkey.x>>i)&0x1,(lkey.y>>i)&0x1,(lkey.z>>i)&0x1,(lkey.w>>i)&0x1 );\n",
        "    flags[index.x+CONFLICT_FREE_OFFSET(index.x)] = 1<<(16*flag.x);\n",
        "    flags[index.y+CONFLICT_FREE_OFFSET(index.y)] = 1<<(16*flag.y);\n",
        "    flags[index.z+CONFLICT_FREE_OFFSET(index.z)] = 1<<(16*flag.z);\n",
        "    flags[index.w+CONFLICT_FREE_OFFSET(index.w)] = 1<<(16*flag.w);\n",
        "    scan (flags);\n",
        "    index.x = (flags[index.x+CONFLICT_FREE_OFFSET(index.x)]>>(16*flag.x))&0xFFFF;\n",
        "    index.y = (flags[index.y+CONFLICT_FREE_OFFSET(index.y)]>>(16*flag.y))&0xFFFF;\n",
        "    index.z = (flags[index.z+CONFLICT_FREE_OFFSET(index.z)]>>(16*flag.z))&0xFFFF;\n",
        "    index.w = (flags[index.w+CONFLICT_FREE_OFFSET(index.w)]>>(16*flag.w))&0xFFFF;\n",
        "    unsigned short offset = flags[4*get_local_size(0)+CONFLICT_FREE_OFFSET(4*get_local_size(0))]&0xFFFF;\n",
        "    index.x += (flag.x) ? offset : 0;\n",
        "    index.y += (flag.y) ? offset : 0;\n",
        "    index.z += (flag.z) ? offset : 0;\n",
        "    index.w += (flag.w) ? offset : 0;\n",
        "    barrier(CLK_LOCAL_MEM_FENCE );\n",
        "  }\n",
        "  if (gid < numElems){\n",
        "    keys[get_group_id(0)*4*SORT_BS+index.x] = lkey.x;\n",
        "    keys[get_group_id(0)*4*SORT_BS+index.y] = lkey.y;\n",
        "    keys[get_group_id(0)*4*SORT_BS+index.z] = lkey.z;\n",
        "    keys[get_group_id(0)*4*SORT_BS+index.w] = lkey.w;\n",
        "    values[get_group_id(0)*4*SORT_BS+index.x] = lvalue.x;\n",
        "    values[get_group_id(0)*4*SORT_BS+index.y] = lvalue.y;\n",
        "    values[get_group_id(0)*4*SORT_BS+index.z] = lvalue.z;\n",
        "    values[get_group_id(0)*4*SORT_BS+index.w] = lvalue.w;\n",
        "  }\n",
        "  if (tid < (1<<BITS)){\n",
        "    histo[get_num_groups(0)*get_local_id(0)+get_group_id(0)] = histo_s[tid];\n",
        "  }\n",
        "}\n",
        "\n",
        "__kernel void splitRearrange (int numElems, int iter,\n",
        "  __global unsigned int* keys_i,\n",
        "  __global unsigned int* keys_o,\n",
        "  __global unsigned int* values_i,\n",
        "  __global unsigned int* values_o,\n",
        "  __global unsigned int* histo){\n",
        "  __local unsigned int histo_s[(1<<BITS)];\n",
        "  __local uint array_s[4*SORT_BS];\n",
        "  int index = get_group_id(0)*4*SORT_BS + 4*get_local_id(0);\n",
        "  if (get_local_id(0) < (1<<BITS)){\n",
        "    histo_s[get_local_id(0)] = histo[get_num_groups(0)*get_local_id(0)+get_group_id(0)];\n",
        "  }\n",
        "  uint4 mine, value;\n",
        "  if (index < numElems){\n",
        "    mine = *((__global uint4*)(keys_i+index));\n",
        "    value = *((__global uint4*)(values_i+index));\n",
        "  } else {\n",
        "    mine.x = UINT32_MAX;\n",
        "    mine.y = UINT32_MAX;\n",
        "    mine.z = UINT32_MAX;\n",
        "    mine.w = UINT32_MAX;\n",
        "  }\n",
        "  uint4 masks = (uint4) ( (mine.x&((1<<(BITS*(iter+1)))-1))>>(BITS*iter),\n",
        "   (mine.y&((1<<(BITS*(iter+1)))-1))>>(BITS*iter),\n",
        "   (mine.z&((1<<(BITS*(iter+1)))-1))>>(BITS*iter),\n",
        "   (mine.w&((1<<(BITS*(iter+1)))-1))>>(BITS*iter) );\n",
        "  vstore4(masks, get_local_id(0), (__local uint *)array_s);\n",
        "  barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "  uint4 new_index = (uint4) ( histo_s[masks.x],histo_s[masks.y],histo_s[masks.z],histo_s[masks.w] );\n",
        "  int i = 4*get_local_id(0)-1;\n",
        "  while (i >= 0){\n",
        "    if (array_s[i] == masks.x){\n",
        "      new_index.x++;\n",
        "      i--;\n",
        "    } else {\n",
        "      break;\n",
        "    }\n",
        "  }\n",
        "  new_index.y = (masks.y == masks.x) ? new_index.x+1 : new_index.y;\n",
        "  new_index.z = (masks.z == masks.y) ? new_index.y+1 : new_index.z;\n",
        "  new_index.w = (masks.w == masks.z) ? new_index.z+1 : new_index.w;\n",
        "  if (index < numElems){\n",
        "    keys_o[new_index.x] = mine.x;\n",
        "    values_o[new_index.x] = value.x;\n",
        "\n",
        "    keys_o[new_index.y] = mine.y;\n",
        "    values_o[new_index.y] = value.y;\n",
        "\n",
        "    keys_o[new_index.z] = mine.z;\n",
        "    values_o[new_index.z] = value.z;\n",
        "\n",
        "    keys_o[new_index.w] = mine.w;\n",
        "    values_o[new_index.w] = value.w;\n",
        "  }\n",
        "}\n",
        "\n",
        "#define BLOCK_SIZE 1024\n",
        "#define GRID_SIZE 65535\n",
        "#define NUM_BANKS 16\n",
        "#define LOG_NUM_BANKS 4\n",
        "#define LNB LOG_NUM_BANKS\n",
        "#define EXPANDED_SIZE(__x) (__x+(__x>>LOG_NUM_BANKS)+(__x>>(2*LOG_NUM_BANKS)))\n",
        "__kernel void scan_L1_kernel(unsigned int n, __global unsigned int* dataBase, unsigned int data_offset, __global unsigned int* interBase, unsigned int inter_offset)\n",
        "{\n",
        "  __local unsigned int s_data[EXPANDED_SIZE(BLOCK_SIZE)];\n",
        "  __global unsigned int *data = dataBase + data_offset;\n",
        "  __global unsigned int *inter = interBase + inter_offset;\n",
        "  unsigned int thid = get_local_id(0);\n",
        "  unsigned int g_ai = get_group_id(0)*2*get_local_size(0) + get_local_id(0);\n",
        "  unsigned int g_bi = g_ai + get_local_size(0);\n",
        "  unsigned int s_ai = thid;\n",
        "  unsigned int s_bi = thid + get_local_size(0);\n",
        "  s_ai += CONFLICT_FREE_OFFSET(s_ai);\n",
        "  s_bi += CONFLICT_FREE_OFFSET(s_bi);\n",
        "  s_data[s_ai] = (g_ai < n) ? data[g_ai] : 0;\n",
        "  s_data[s_bi] = (g_bi < n) ? data[g_bi] : 0;\n",
        "  unsigned int stride = 1;\n",
        "  for (unsigned int d = get_local_size(0); d > 0; d >>= 1) {\n",
        "      barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "      if (thid < d) {\n",
        "        unsigned int i  = 2*stride*thid;\n",
        "        unsigned int ai = i + stride - 1;\n",
        "        unsigned int bi = ai + stride;\n",
        "        ai += CONFLICT_FREE_OFFSET(ai);\n",
        "        bi += CONFLICT_FREE_OFFSET(bi);\n",
        "        s_data[bi] += s_data[ai];\n",
        "      }\n",
        "      stride *= 2;\n",
        "    }\n",
        "    if (thid == 0) {\n",
        "      unsigned int last = get_local_size(0)*2 -1;\n",
        "      last += CONFLICT_FREE_OFFSET(last);\n",
        "      inter[get_group_id(0)] = s_data[last];\n",
        "      s_data[last] = 0;\n",
        "    }\n",
        "    for (unsigned int d = 1; d <= get_local_size(0); d *= 2) {\n",
        "      stride >>= 1;\n",
        "      barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "      if (thid < d) {\n",
        "        unsigned int i  = 2*stride*thid;\n",
        "        unsigned int ai = i + stride - 1;\n",
        "        unsigned int bi = ai + stride;\n",
        "\n",
        "        ai += CONFLICT_FREE_OFFSET(ai);\n",
        "        bi += CONFLICT_FREE_OFFSET(bi);\n",
        "\n",
        "        unsigned int t  = s_data[ai];\n",
        "        s_data[ai] = s_data[bi];\n",
        "        s_data[bi] += t;\n",
        "      }\n",
        "    }\n",
        "    barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "    if (g_ai < n) { data[g_ai] = s_data[s_ai]; }\n",
        "    if (g_bi < n) { data[g_bi] = s_data[s_bi]; }\n",
        "  }\n",
        "\n",
        "  __kernel void scan_inter1_kernel(__global unsigned int* data, unsigned int iter)\n",
        "  {\n",
        "    __local unsigned int s_data[DYN_LOCAL_MEM_SIZE];\n",
        "    unsigned int thid = get_local_id(0);\n",
        "    unsigned int gthid = get_global_id(0);\n",
        "    unsigned int gi = 2*iter*gthid;\n",
        "    unsigned int g_ai = gi + iter - 1;\n",
        "    unsigned int g_bi = g_ai + iter;\n",
        "    unsigned int s_ai = 2*thid;\n",
        "    unsigned int s_bi = 2*thid + 1;\n",
        "    s_ai += CONFLICT_FREE_OFFSET(s_ai);\n",
        "    s_bi += CONFLICT_FREE_OFFSET(s_bi);\n",
        "    s_data[s_ai] = data[g_ai];\n",
        "    s_data[s_bi] = data[g_bi];\n",
        "    unsigned int stride = 1;\n",
        "    for (unsigned int d = get_local_size(0); d > 0; d >>= 1) {\n",
        "      barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "      if (thid < d) {\n",
        "        unsigned int i  = 2*stride*thid;\n",
        "        unsigned int ai = i + stride - 1;\n",
        "        unsigned int bi = ai + stride;\n",
        "\n",
        "        ai += CONFLICT_FREE_OFFSET(ai);\n",
        "        bi += CONFLICT_FREE_OFFSET(bi);\n",
        "        s_data[bi] += s_data[ai];\n",
        "      }\n",
        "      stride *= 2;\n",
        "    }\n",
        "    barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "    data[g_ai] = s_data[s_ai];\n",
        "    data[g_bi] = s_data[s_bi];\n",
        "  }\n",
        "\n",
        "  __kernel void scan_inter2_kernel(__global unsigned int* data, unsigned int iter)\n",
        "  {\n",
        "    __local unsigned int s_data[DYN_LOCAL_MEM_SIZE];\n",
        "\n",
        "    unsigned int thid = get_local_id(0);\n",
        "    unsigned int gthid = get_global_id(0);\n",
        "    unsigned int gi = 2*iter*gthid;\n",
        "    unsigned int g_ai = gi + iter - 1;\n",
        "    unsigned int g_bi = g_ai + iter;\n",
        "    unsigned int s_ai = 2*thid;\n",
        "    unsigned int s_bi = 2*thid + 1;\n",
        "    s_ai += CONFLICT_FREE_OFFSET(s_ai);\n",
        "    s_bi += CONFLICT_FREE_OFFSET(s_bi);\n",
        "    s_data[s_ai] = data[g_ai];\n",
        "    s_data[s_bi] = data[g_bi];\n",
        "    unsigned int stride = get_local_size(0)*2;\n",
        "    for (unsigned int d = 1; d <= get_local_size(0); d *= 2) {\n",
        "      stride >>= 1;\n",
        "      barrier(CLK_LOCAL_MEM_FENCE );\n",
        "      if (thid < d) {\n",
        "        unsigned int i  = 2*stride*thid;\n",
        "        unsigned int ai = i + stride - 1;\n",
        "        unsigned int bi = ai + stride;\n",
        "        ai += CONFLICT_FREE_OFFSET(ai);\n",
        "        bi += CONFLICT_FREE_OFFSET(bi);\n",
        "        unsigned int t  = s_data[ai];\n",
        "        s_data[ai] = s_data[bi];\n",
        "        s_data[bi] += t;\n",
        "      }\n",
        "    }\n",
        "    barrier(CLK_LOCAL_MEM_FENCE );\n",
        "    data[g_ai] = s_data[s_ai];\n",
        "    data[g_bi] = s_data[s_bi];\n",
        "  }\n",
        "\n",
        "__kernel void uniformAdd(unsigned int n, __global unsigned int *dataBase, unsigned int data_offset, __global unsigned int *interBase, unsigned int inter_offset){\n",
        "    __local unsigned int uni;\n",
        "    __global unsigned int *data = dataBase + data_offset;\n",
        "    __global unsigned int *inter = interBase + inter_offset;\n",
        "    if (get_local_id(0) == 0) { uni = inter[get_group_id(0)]; }\n",
        "    barrier(CLK_LOCAL_MEM_FENCE ); //__syncthreads();\n",
        "    unsigned int g_ai = get_group_id(0)*2*get_local_size(0) + get_local_id(0);\n",
        "    unsigned int g_bi = g_ai + get_local_size(0);\n",
        "    if (g_ai < n) { data[g_ai] += uni; }\n",
        "    if (g_bi < n) { data[g_bi] += uni; }\n",
        "}\n",
        "\n",
        "\n"
};


#endif //MOBILEHETEROGENOUSPROJECT_CLION_MRI_GRIDDING_H
