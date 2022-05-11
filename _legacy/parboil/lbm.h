//
// Created by pfxu on 20/8/2019.
//

#ifndef MOBILEHETEROGENOUSPROJECT_CLION_LBM_H
#define MOBILEHETEROGENOUSPROJECT_CLION_LBM_H

#include "benchmark.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <sys/stat.h>

typedef enum {
    OBSTACLE = 1 << 0,
    ACCEL = 1 << 1,
    IN_OUT_FLOW = 1 << 2
} CELL_FLAGS;

#define N_TIME_STEPS 32

#define OMEGA (1.95f)

#define OUTPUT_PRECISION float

typedef enum {
    C = 0,
    N, S, E, W, T, B,
    NE, NW, SE, SW,
    NT, NB, ST, SB,
    ET, EB, WT, WB,
    FLAGS, N_CELL_ENTRIES
} CELL_ENTRIES;

typedef enum {
    NOTHING = 0, COMPARE, STORE
} MAIN_Action;
typedef enum {
    LDC = 0, CHANNEL
} MAIN_SimType;

typedef struct {
    int nTimeSteps;
    char const *resultFilename;
    MAIN_Action action;
    MAIN_SimType simType;
    char const *obstacleFilename;
} MAIN_Param;

#define SIZE   (120)
#define SIZE_X (1*SIZE)
#define SIZE_Y (1*SIZE)
#define SIZE_Z (150)
/*############################################################################*/

typedef float LBM_Grid[SIZE_Z * SIZE_Y * SIZE_X * N_CELL_ENTRIES];
typedef LBM_Grid *LBM_GridPtr;

typedef float *GLBM_Grid;//float LBM_Grid[PADDED_Z*PADDED_Y*PADDED_X*N_CELL_ENTRIES];
typedef GLBM_Grid *GLBM_GridPtr;

/*############################################################################*/

#define CALC_INDEX(x, y, z, e) ((e)+N_CELL_ENTRIES*((x)+ \
                         (y)*SIZE_X+(z)*SIZE_X*SIZE_Y))

#define SWEEP_VAR int i;

#define SWEEP_START(x1, y1, z1, x2, y2, z2) \
  for( i = CALC_INDEX(x1, y1, z1, 0); \
       i < CALC_INDEX(x2, y2, z2, 0); \
       i += N_CELL_ENTRIES ) {

#define SWEEP_END }

#define SWEEP_X  ((i / N_CELL_ENTRIES) % SIZE_X)
#define SWEEP_Y (((i / N_CELL_ENTRIES) / SIZE_X) % SIZE_Y)
#define SWEEP_Z  ((i / N_CELL_ENTRIES) / (SIZE_X*SIZE_Y))
#define GRID_ENTRY(g, x, y, z, e)          ((g)[CALC_INDEX( x,  y,  z, e)])
#define GRID_ENTRY_SWEEP(g, dx, dy, dz, e) ((g)[CALC_INDEX(dx, dy, dz, e)+(i)])
#define LOCAL(g, e)       (GRID_ENTRY_SWEEP( g,  0,  0,  0, e ))
#define NEIGHBOR_C(g, e)  (GRID_ENTRY_SWEEP( g,  0,  0,  0, e ))
#define NEIGHBOR_N(g, e)  (GRID_ENTRY_SWEEP( g,  0, +1,  0, e ))
#define NEIGHBOR_S(g, e)  (GRID_ENTRY_SWEEP( g,  0, -1,  0, e ))
#define NEIGHBOR_E(g, e)  (GRID_ENTRY_SWEEP( g, +1,  0,  0, e ))
#define NEIGHBOR_W(g, e)  (GRID_ENTRY_SWEEP( g, -1,  0,  0, e ))
#define NEIGHBOR_T(g, e)  (GRID_ENTRY_SWEEP( g,  0,  0, +1, e ))
#define NEIGHBOR_B(g, e)  (GRID_ENTRY_SWEEP( g,  0,  0, -1, e ))
#define NEIGHBOR_NE(g, e) (GRID_ENTRY_SWEEP( g, +1, +1,  0, e ))
#define NEIGHBOR_NW(g, e) (GRID_ENTRY_SWEEP( g, -1, +1,  0, e ))
#define NEIGHBOR_SE(g, e) (GRID_ENTRY_SWEEP( g, +1, -1,  0, e ))
#define NEIGHBOR_SW(g, e) (GRID_ENTRY_SWEEP( g, -1, -1,  0, e ))
#define NEIGHBOR_NT(g, e) (GRID_ENTRY_SWEEP( g,  0, +1, +1, e ))
#define NEIGHBOR_NB(g, e) (GRID_ENTRY_SWEEP( g,  0, +1, -1, e ))
#define NEIGHBOR_ST(g, e) (GRID_ENTRY_SWEEP( g,  0, -1, +1, e ))
#define NEIGHBOR_SB(g, e) (GRID_ENTRY_SWEEP( g,  0, -1, -1, e ))
#define NEIGHBOR_ET(g, e) (GRID_ENTRY_SWEEP( g, +1,  0, +1, e ))
#define NEIGHBOR_EB(g, e) (GRID_ENTRY_SWEEP( g, +1,  0, -1, e ))
#define NEIGHBOR_WT(g, e) (GRID_ENTRY_SWEEP( g, -1,  0, +1, e ))
#define NEIGHBOR_WB(g, e) (GRID_ENTRY_SWEEP( g, -1,  0, -1, e ))
#define COLLIDE_STREAM
#ifdef COLLIDE_STREAM
#define SRC_C(g)  (LOCAL( g, C  ))
#define SRC_N(g)  (LOCAL( g, N  ))
#define SRC_S(g)  (LOCAL( g, S  ))
#define SRC_E(g)  (LOCAL( g, E  ))
#define SRC_W(g)  (LOCAL( g, W  ))
#define SRC_T(g)  (LOCAL( g, T  ))
#define SRC_B(g)  (LOCAL( g, B  ))
#define SRC_NE(g) (LOCAL( g, NE ))
#define SRC_NW(g) (LOCAL( g, NW ))
#define SRC_SE(g) (LOCAL( g, SE ))
#define SRC_SW(g) (LOCAL( g, SW ))
#define SRC_NT(g) (LOCAL( g, NT ))
#define SRC_NB(g) (LOCAL( g, NB ))
#define SRC_ST(g) (LOCAL( g, ST ))
#define SRC_SB(g) (LOCAL( g, SB ))
#define SRC_ET(g) (LOCAL( g, ET ))
#define SRC_EB(g) (LOCAL( g, EB ))
#define SRC_WT(g) (LOCAL( g, WT ))
#define SRC_WB(g) (LOCAL( g, WB ))
#define DST_C(g)  (NEIGHBOR_C ( g, C  ))
#define DST_N(g)  (NEIGHBOR_N ( g, N  ))
#define DST_S(g)  (NEIGHBOR_S ( g, S  ))
#define DST_E(g)  (NEIGHBOR_E ( g, E  ))
#define DST_W(g)  (NEIGHBOR_W ( g, W  ))
#define DST_T(g)  (NEIGHBOR_T ( g, T  ))
#define DST_B(g)  (NEIGHBOR_B ( g, B  ))
#define DST_NE(g) (NEIGHBOR_NE( g, NE ))
#define DST_NW(g) (NEIGHBOR_NW( g, NW ))
#define DST_SE(g) (NEIGHBOR_SE( g, SE ))
#define DST_SW(g) (NEIGHBOR_SW( g, SW ))
#define DST_NT(g) (NEIGHBOR_NT( g, NT ))
#define DST_NB(g) (NEIGHBOR_NB( g, NB ))
#define DST_ST(g) (NEIGHBOR_ST( g, ST ))
#define DST_SB(g) (NEIGHBOR_SB( g, SB ))
#define DST_ET(g) (NEIGHBOR_ET( g, ET ))
#define DST_EB(g) (NEIGHBOR_EB( g, EB ))
#define DST_WT(g) (NEIGHBOR_WT( g, WT ))
#define DST_WB(g) (NEIGHBOR_WB( g, WB ))
#else /* COLLIDE_STREAM */
#define SRC_C(g)  (NEIGHBOR_C ( g, C  ))
#define SRC_N(g)  (NEIGHBOR_S ( g, N  ))
#define SRC_S(g)  (NEIGHBOR_N ( g, S  ))
#define SRC_E(g)  (NEIGHBOR_W ( g, E  ))
#define SRC_W(g)  (NEIGHBOR_E ( g, W  ))
#define SRC_T(g)  (NEIGHBOR_B ( g, T  ))
#define SRC_B(g)  (NEIGHBOR_T ( g, B  ))
#define SRC_NE(g) (NEIGHBOR_SW( g, NE ))
#define SRC_NW(g) (NEIGHBOR_SE( g, NW ))
#define SRC_SE(g) (NEIGHBOR_NW( g, SE ))
#define SRC_SW(g) (NEIGHBOR_NE( g, SW ))
#define SRC_NT(g) (NEIGHBOR_SB( g, NT ))
#define SRC_NB(g) (NEIGHBOR_ST( g, NB ))
#define SRC_ST(g) (NEIGHBOR_NB( g, ST ))
#define SRC_SB(g) (NEIGHBOR_NT( g, SB ))
#define SRC_ET(g) (NEIGHBOR_WB( g, ET ))
#define SRC_EB(g) (NEIGHBOR_WT( g, EB ))
#define SRC_WT(g) (NEIGHBOR_EB( g, WT ))
#define SRC_WB(g) (NEIGHBOR_ET( g, WB ))
#define DST_C(g)  (LOCAL( g, C  ))
#define DST_N(g)  (LOCAL( g, N  ))
#define DST_S(g)  (LOCAL( g, S  ))
#define DST_E(g)  (LOCAL( g, E  ))
#define DST_W(g)  (LOCAL( g, W  ))
#define DST_T(g)  (LOCAL( g, T  ))
#define DST_B(g)  (LOCAL( g, B  ))
#define DST_NE(g) (LOCAL( g, NE ))
#define DST_NW(g) (LOCAL( g, NW ))
#define DST_SE(g) (LOCAL( g, SE ))
#define DST_SW(g) (LOCAL( g, SW ))
#define DST_NT(g) (LOCAL( g, NT ))
#define DST_NB(g) (LOCAL( g, NB ))
#define DST_ST(g) (LOCAL( g, ST ))
#define DST_SB(g) (LOCAL( g, SB ))
#define DST_ET(g) (LOCAL( g, ET ))
#define DST_EB(g) (LOCAL( g, EB ))
#define DST_WT(g) (LOCAL( g, WT ))
#define DST_WB(g) (LOCAL( g, WB ))
#endif /* COLLIDE_STREAM */

#define MAGIC_CAST(v) ((unsigned int*) ((void*) (&(v))))
#define FLAG_VAR(v) unsigned int* const _aux_ = MAGIC_CAST(v)
#define TEST_FLAG_SWEEP(g, f)     ((*MAGIC_CAST(LOCAL(g, FLAGS))) & (f))
#define SET_FLAG_SWEEP(g, f)      {FLAG_VAR(LOCAL(g, FLAGS)); (*_aux_) |=  (f);}
#define CLEAR_FLAG_SWEEP(g, f)    {FLAG_VAR(LOCAL(g, FLAGS)); (*_aux_) &= ~(f);}
#define CLEAR_ALL_FLAGS_SWEEP(g) {FLAG_VAR(LOCAL(g, FLAGS)); (*_aux_)  =    0;}
#define TEST_FLAG(g, x, y, z, f)     ((*MAGIC_CAST(GRID_ENTRY(g, x, y, z, FLAGS))) & (f))
#define SET_FLAG(g, x, y, z, f)      {FLAG_VAR(GRID_ENTRY(g, x, y, z, FLAGS)); (*_aux_) |=  (f);}
#define CLEAR_FLAG(g, x, y, z, f)    {FLAG_VAR(GRID_ENTRY(g, x, y, z, FLAGS)); (*_aux_) &= ~(f);}
#define CLEAR_ALL_FLAGS(g, x, y, z) {FLAG_VAR(GRID_ENTRY(g, x, y, z, FLAGS)); (*_aux_)  =    0;}

#define DFL1 (1.0/ 3.0)
#define DFL2 (1.0/18.0)
#define DFL3 (1.0/36.0)

/*############################################################################*/

void LBM_allocateGrid(float **ptr) {
    const size_t margin = 2 * SIZE_X * SIZE_Y * N_CELL_ENTRIES, size = sizeof(LBM_Grid) + 2 * margin * sizeof(float);
    *ptr = static_cast<float *>(malloc(size));
    if (!*ptr) {
        printf("LBM_allocateGrid: could not allocate %.1f MByte\n",
               size / (1024.0 * 1024.0));
        exit(1);
    }
#if !defined(SPEC_CPU)
    printf("LBM_allocateGrid: allocated %.1f MByte\n", size / (1024.0 * 1024.0));
#endif
    *ptr += margin;
}

/*############################################################################*/

void LBM_freeGrid(float **ptr) {
    const size_t margin = 2 * SIZE_X * SIZE_Y * N_CELL_ENTRIES;
    free(*ptr - margin);
    *ptr = NULL;
}

/*############################################################################*/

void LBM_initializeGrid(LBM_Grid grid) {
    SWEEP_VAR
#pragma omp parallel for
    SWEEP_START(0, 0, -2, 0, 0, SIZE_Z + 2)
        LOCAL(grid, C) = DFL1;
        LOCAL(grid, N) = DFL2;
        LOCAL(grid, S) = DFL2;
        LOCAL(grid, E) = DFL2;
        LOCAL(grid, W) = DFL2;
        LOCAL(grid, T) = DFL2;
        LOCAL(grid, B) = DFL2;
        LOCAL(grid, NE) = DFL3;
        LOCAL(grid, NW) = DFL3;
        LOCAL(grid, SE) = DFL3;
        LOCAL(grid, SW) = DFL3;
        LOCAL(grid, NT) = DFL3;
        LOCAL(grid, NB) = DFL3;
        LOCAL(grid, ST) = DFL3;
        LOCAL(grid, SB) = DFL3;
        LOCAL(grid, ET) = DFL3;
        LOCAL(grid, EB) = DFL3;
        LOCAL(grid, WT) = DFL3;
        LOCAL(grid, WB) = DFL3;

        CLEAR_ALL_FLAGS_SWEEP(grid);
    SWEEP_END
}

/*############################################################################*/

void LBM_swapGrids(LBM_GridPtr *grid1, LBM_GridPtr *grid2) {
    LBM_GridPtr aux = *grid1;
    *grid1 = *grid2;
    *grid2 = aux;
}

void LBM_loadObstacleFile(LBM_Grid grid, const char *filename) {
    int x, y, z;
    FILE *file = fopen(filename, "rb");
    for (z = 0; z < SIZE_Z; z++) {
        for (y = 0; y < SIZE_Y; y++) {
            for (x = 0; x < SIZE_X; x++) {
                if (fgetc(file) != '.') SET_FLAG(grid, x, y, z, OBSTACLE);
            }
            fgetc(file);
        }
        fgetc(file);
    }
    fclose(file);
}

/*############################################################################*/

void LBM_initializeSpecialCellsForLDC(LBM_Grid grid) {
    int x, y, z;
#pragma omp parallel for private( x, y )
    for (z = -2; z < SIZE_Z + 2; z++) {
        for (y = 0; y < SIZE_Y; y++) {
            for (x = 0; x < SIZE_X; x++) {
                if (x == 0 || x == SIZE_X - 1 ||
                    y == 0 || y == SIZE_Y - 1 ||
                    z == 0 || z == SIZE_Z - 1) {
                    SET_FLAG(grid, x, y, z, OBSTACLE);
                } else {
                    if ((z == 1 || z == SIZE_Z - 2) &&
                        x > 1 && x < SIZE_X - 2 &&
                        y > 1 && y < SIZE_Y - 2) {
                        SET_FLAG(grid, x, y, z, ACCEL);
                    }
                }
            }
        }
    }
}

void LBM_initializeSpecialCellsForChannel(LBM_Grid grid) {
    int x, y, z;
#pragma omp parallel for private( x, y )
    for (z = -2; z < SIZE_Z + 2; z++) {
        for (y = 0; y < SIZE_Y; y++) {
            for (x = 0; x < SIZE_X; x++) {
                if (x == 0 || x == SIZE_X - 1 ||
                    y == 0 || y == SIZE_Y - 1) {
                    SET_FLAG(grid, x, y, z, OBSTACLE);

                    if ((z == 0 || z == SIZE_Z - 1) &&
                        !TEST_FLAG(grid, x, y, z, OBSTACLE)) SET_FLAG(grid, x, y, z, IN_OUT_FLOW);
                }
            }
        }
    }
}

void LBM_performStreamCollide(LBM_Grid srcGrid, LBM_Grid dstGrid) {
    SWEEP_VAR
    float ux, uy, uz, u2, rho;
#pragma omp parallel for private( ux, uy, uz, u2, rho )
    SWEEP_START(0, 0, 0, 0, 0, SIZE_Z)
        if (TEST_FLAG_SWEEP(srcGrid, OBSTACLE)) {
            DST_C (dstGrid) = SRC_C (srcGrid);
            DST_S (dstGrid) = SRC_N (srcGrid);
            DST_N (dstGrid) = SRC_S (srcGrid);
            DST_W (dstGrid) = SRC_E (srcGrid);
            DST_E (dstGrid) = SRC_W (srcGrid);
            DST_B (dstGrid) = SRC_T (srcGrid);
            DST_T (dstGrid) = SRC_B (srcGrid);
            DST_SW(dstGrid) = SRC_NE(srcGrid);
            DST_SE(dstGrid) = SRC_NW(srcGrid);
            DST_NW(dstGrid) = SRC_SE(srcGrid);
            DST_NE(dstGrid) = SRC_SW(srcGrid);
            DST_SB(dstGrid) = SRC_NT(srcGrid);
            DST_ST(dstGrid) = SRC_NB(srcGrid);
            DST_NB(dstGrid) = SRC_ST(srcGrid);
            DST_NT(dstGrid) = SRC_SB(srcGrid);
            DST_WB(dstGrid) = SRC_ET(srcGrid);
            DST_WT(dstGrid) = SRC_EB(srcGrid);
            DST_EB(dstGrid) = SRC_WT(srcGrid);
            DST_ET(dstGrid) = SRC_WB(srcGrid);
            continue;
        }
        rho = +SRC_C (srcGrid) + SRC_N (srcGrid)
              + SRC_S (srcGrid) + SRC_E (srcGrid)
              + SRC_W (srcGrid) + SRC_T (srcGrid)
              + SRC_B (srcGrid) + SRC_NE(srcGrid)
              + SRC_NW(srcGrid) + SRC_SE(srcGrid)
              + SRC_SW(srcGrid) + SRC_NT(srcGrid)
              + SRC_NB(srcGrid) + SRC_ST(srcGrid)
              + SRC_SB(srcGrid) + SRC_ET(srcGrid)
              + SRC_EB(srcGrid) + SRC_WT(srcGrid)
              + SRC_WB(srcGrid);
        ux = +SRC_E (srcGrid) - SRC_W (srcGrid)
             + SRC_NE(srcGrid) - SRC_NW(srcGrid)
             + SRC_SE(srcGrid) - SRC_SW(srcGrid)
             + SRC_ET(srcGrid) + SRC_EB(srcGrid)
             - SRC_WT(srcGrid) - SRC_WB(srcGrid);
        uy = +SRC_N (srcGrid) - SRC_S (srcGrid)
             + SRC_NE(srcGrid) + SRC_NW(srcGrid)
             - SRC_SE(srcGrid) - SRC_SW(srcGrid)
             + SRC_NT(srcGrid) + SRC_NB(srcGrid)
             - SRC_ST(srcGrid) - SRC_SB(srcGrid);
        uz = +SRC_T (srcGrid) - SRC_B (srcGrid)
             + SRC_NT(srcGrid) - SRC_NB(srcGrid)
             + SRC_ST(srcGrid) - SRC_SB(srcGrid)
             + SRC_ET(srcGrid) - SRC_EB(srcGrid)
             + SRC_WT(srcGrid) - SRC_WB(srcGrid);
        ux /= rho;
        uy /= rho;
        uz /= rho;
        if (TEST_FLAG_SWEEP(srcGrid, ACCEL)) {
            ux = 0.005f;
            uy = 0.002f;
            uz = 0.000f;
        }
        u2 = 1.5f * (ux * ux + uy * uy + uz * uz);
        DST_C (dstGrid) = (1.0f - OMEGA) * SRC_C (srcGrid) + DFL1 * OMEGA * rho * (1.0f - u2);
        DST_N (dstGrid) = (1.0f - OMEGA) * SRC_N (srcGrid) + DFL2 * OMEGA * rho * (1.0f + uy * (4.5f * uy + 3.0f) - u2);
        DST_S (dstGrid) = (1.0f - OMEGA) * SRC_S (srcGrid) + DFL2 * OMEGA * rho * (1.0f + uy * (4.5f * uy - 3.0f) - u2);
        DST_E (dstGrid) = (1.0f - OMEGA) * SRC_E (srcGrid) + DFL2 * OMEGA * rho * (1.0f + ux * (4.5f * ux + 3.0f) - u2);
        DST_W (dstGrid) = (1.0f - OMEGA) * SRC_W (srcGrid) + DFL2 * OMEGA * rho * (1.0f + ux * (4.5f * ux - 3.0f) - u2);
        DST_T (dstGrid) = (1.0f - OMEGA) * SRC_T (srcGrid) + DFL2 * OMEGA * rho * (1.0f + uz * (4.5f * uz + 3.0f) - u2);
        DST_B (dstGrid) = (1.0f - OMEGA) * SRC_B (srcGrid) + DFL2 * OMEGA * rho * (1.0f + uz * (4.5f * uz - 3.0f) - u2);
        DST_NE(dstGrid) = (1.0f - OMEGA) * SRC_NE(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (+ux + uy) * (4.5f * (+ux + uy) + 3.0f) - u2);
        DST_NW(dstGrid) = (1.0f - OMEGA) * SRC_NW(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (-ux + uy) * (4.5f * (-ux + uy) + 3.0f) - u2);
        DST_SE(dstGrid) = (1.0f - OMEGA) * SRC_SE(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (+ux - uy) * (4.5f * (+ux - uy) + 3.0f) - u2);
        DST_SW(dstGrid) = (1.0f - OMEGA) * SRC_SW(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (-ux - uy) * (4.5f * (-ux - uy) + 3.0f) - u2);
        DST_NT(dstGrid) = (1.0f - OMEGA) * SRC_NT(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (+uy + uz) * (4.5f * (+uy + uz) + 3.0f) - u2);
        DST_NB(dstGrid) = (1.0f - OMEGA) * SRC_NB(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (+uy - uz) * (4.5f * (+uy - uz) + 3.0f) - u2);
        DST_ST(dstGrid) = (1.0f - OMEGA) * SRC_ST(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (-uy + uz) * (4.5f * (-uy + uz) + 3.0f) - u2);
        DST_SB(dstGrid) = (1.0f - OMEGA) * SRC_SB(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (-uy - uz) * (4.5f * (-uy - uz) + 3.0f) - u2);
        DST_ET(dstGrid) = (1.0f - OMEGA) * SRC_ET(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (+ux + uz) * (4.5f * (+ux + uz) + 3.0f) - u2);
        DST_EB(dstGrid) = (1.0f - OMEGA) * SRC_EB(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (+ux - uz) * (4.5f * (+ux - uz) + 3.0f) - u2);
        DST_WT(dstGrid) = (1.0f - OMEGA) * SRC_WT(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (-ux + uz) * (4.5f * (-ux + uz) + 3.0f) - u2);
        DST_WB(dstGrid) = (1.0f - OMEGA) * SRC_WB(srcGrid) +
                          DFL3 * OMEGA * rho * (1.0f + (-ux - uz) * (4.5f * (-ux - uz) + 3.0f) - u2);
    SWEEP_END
}

void LBM_handleInOutFlow(LBM_Grid srcGrid) {
    float ux, uy, uz, rho,
            ux1, uy1, uz1, rho1,
            ux2, uy2, uz2, rho2,
            u2, px, py;
    SWEEP_VAR

#pragma omp parallel for private( ux, uy, uz, rho, ux1, uy1, uz1, rho1, ux2, uy2, uz2, rho2, u2, px, py )
    SWEEP_START(0, 0, 0, 0, 0, 1)
        rho1 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, C) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, N)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, S) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, E)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, W) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, T)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, B) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, NE)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, NW) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, SE)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, SW) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, NT)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, NB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, ST)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, SB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, ET)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, EB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, WT)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 1, WB);
        rho2 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, C) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, N)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, S) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, E)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, W) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, T)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, B) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, NE)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, NW) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, SE)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, SW) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, NT)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, NB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, ST)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, SB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, ET)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, EB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, WT)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, 2, WB);
        rho = 2.0 * rho1 - rho2;
        px = (SWEEP_X / (0.5 * (SIZE_X - 1))) - 1.0;
        py = (SWEEP_Y / (0.5 * (SIZE_Y - 1))) - 1.0;
        ux = 0.00;
        uy = 0.00;
        uz = 0.01 * (1.0 - px * px) * (1.0 - py * py);
        u2 = 1.5 * (ux * ux + uy * uy + uz * uz);
        LOCAL(srcGrid, C) = DFL1 * rho * (1.0 - u2);
        LOCAL(srcGrid, N) = DFL2 * rho * (1.0 + uy * (4.5 * uy + 3.0) - u2);
        LOCAL(srcGrid, S) = DFL2 * rho * (1.0 + uy * (4.5 * uy - 3.0) - u2);
        LOCAL(srcGrid, E) = DFL2 * rho * (1.0 + ux * (4.5 * ux + 3.0) - u2);
        LOCAL(srcGrid, W) = DFL2 * rho * (1.0 + ux * (4.5 * ux - 3.0) - u2);
        LOCAL(srcGrid, T) = DFL2 * rho * (1.0 + uz * (4.5 * uz + 3.0) - u2);
        LOCAL(srcGrid, B) = DFL2 * rho * (1.0 + uz * (4.5 * uz - 3.0) - u2);
        LOCAL(srcGrid, NE) = DFL3 * rho * (1.0 + (+ux + uy) * (4.5 * (+ux + uy) + 3.0) - u2);
        LOCAL(srcGrid, NW) = DFL3 * rho * (1.0 + (-ux + uy) * (4.5 * (-ux + uy) + 3.0) - u2);
        LOCAL(srcGrid, SE) = DFL3 * rho * (1.0 + (+ux - uy) * (4.5 * (+ux - uy) + 3.0) - u2);
        LOCAL(srcGrid, SW) = DFL3 * rho * (1.0 + (-ux - uy) * (4.5 * (-ux - uy) + 3.0) - u2);
        LOCAL(srcGrid, NT) = DFL3 * rho * (1.0 + (+uy + uz) * (4.5 * (+uy + uz) + 3.0) - u2);
        LOCAL(srcGrid, NB) = DFL3 * rho * (1.0 + (+uy - uz) * (4.5 * (+uy - uz) + 3.0) - u2);
        LOCAL(srcGrid, ST) = DFL3 * rho * (1.0 + (-uy + uz) * (4.5 * (-uy + uz) + 3.0) - u2);
        LOCAL(srcGrid, SB) = DFL3 * rho * (1.0 + (-uy - uz) * (4.5 * (-uy - uz) + 3.0) - u2);
        LOCAL(srcGrid, ET) = DFL3 * rho * (1.0 + (+ux + uz) * (4.5 * (+ux + uz) + 3.0) - u2);
        LOCAL(srcGrid, EB) = DFL3 * rho * (1.0 + (+ux - uz) * (4.5 * (+ux - uz) + 3.0) - u2);
        LOCAL(srcGrid, WT) = DFL3 * rho * (1.0 + (-ux + uz) * (4.5 * (-ux + uz) + 3.0) - u2);
        LOCAL(srcGrid, WB) = DFL3 * rho * (1.0 + (-ux - uz) * (4.5 * (-ux - uz) + 3.0) - u2);
    SWEEP_END

#pragma omp parallel for private( ux, uy, uz, rho, ux1, uy1, uz1, rho1, ux2, uy2, uz2, rho2, u2, px, py )
    SWEEP_START(0, 0, SIZE_Z - 1, 0, 0, SIZE_Z)
        rho1 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, C) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, N)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, S) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, E)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, W) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, T)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, B) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NE)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NW) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, SE)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, SW) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NT)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, ST)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, SB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, ET)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, EB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, WT)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, WB);
        ux1 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, E) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, W)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NE) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NW)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, SE) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, SW)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, ET) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, EB)
              - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, WT) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, WB);
        uy1 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, N) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, S)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NE) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NW)
              - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, SE) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, SW)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NT) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NB)
              - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, ST) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, SB);
        uz1 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, T) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, B)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NT) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, NB)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, ST) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, SB)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, ET) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, EB)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, WT) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -1, WB);
        ux1 /= rho1;
        uy1 /= rho1;
        uz1 /= rho1;
        rho2 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, C) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, N)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, S) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, E)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, W) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, T)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, B) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NE)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NW) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, SE)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, SW) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NT)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, ST)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, SB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, ET)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, EB) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, WT)
               + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, WB);
        ux2 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, E) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, W)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NE) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NW)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, SE) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, SW)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, ET) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, EB)
              - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, WT) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, WB);
        uy2 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, N) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, S)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NE) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NW)
              - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, SE) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, SW)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NT) + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NB)
              - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, ST) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, SB);
        uz2 = +GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, T) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, B)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NT) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, NB)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, ST) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, SB)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, ET) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, EB)
              + GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, WT) - GRID_ENTRY_SWEEP(srcGrid, 0, 0, -2, WB);
        ux2 /= rho2;
        uy2 /= rho2;
        uz2 /= rho2;
        rho = 1.0;
        ux = 2 * ux1 - ux2;
        uy = 2 * uy1 - uy2;
        uz = 2 * uz1 - uz2;
        u2 = 1.5 * (ux * ux + uy * uy + uz * uz);
        LOCAL(srcGrid, C) = DFL1 * rho * (1.0 - u2);
        LOCAL(srcGrid, N) = DFL2 * rho * (1.0 + uy * (4.5 * uy + 3.0) - u2);
        LOCAL(srcGrid, S) = DFL2 * rho * (1.0 + uy * (4.5 * uy - 3.0) - u2);
        LOCAL(srcGrid, E) = DFL2 * rho * (1.0 + ux * (4.5 * ux + 3.0) - u2);
        LOCAL(srcGrid, W) = DFL2 * rho * (1.0 + ux * (4.5 * ux - 3.0) - u2);
        LOCAL(srcGrid, T) = DFL2 * rho * (1.0 + uz * (4.5 * uz + 3.0) - u2);
        LOCAL(srcGrid, B) = DFL2 * rho * (1.0 + uz * (4.5 * uz - 3.0) - u2);
        LOCAL(srcGrid, NE) = DFL3 * rho * (1.0 + (+ux + uy) * (4.5 * (+ux + uy) + 3.0) - u2);
        LOCAL(srcGrid, NW) = DFL3 * rho * (1.0 + (-ux + uy) * (4.5 * (-ux + uy) + 3.0) - u2);
        LOCAL(srcGrid, SE) = DFL3 * rho * (1.0 + (+ux - uy) * (4.5 * (+ux - uy) + 3.0) - u2);
        LOCAL(srcGrid, SW) = DFL3 * rho * (1.0 + (-ux - uy) * (4.5 * (-ux - uy) + 3.0) - u2);
        LOCAL(srcGrid, NT) = DFL3 * rho * (1.0 + (+uy + uz) * (4.5 * (+uy + uz) + 3.0) - u2);
        LOCAL(srcGrid, NB) = DFL3 * rho * (1.0 + (+uy - uz) * (4.5 * (+uy - uz) + 3.0) - u2);
        LOCAL(srcGrid, ST) = DFL3 * rho * (1.0 + (-uy + uz) * (4.5 * (-uy + uz) + 3.0) - u2);
        LOCAL(srcGrid, SB) = DFL3 * rho * (1.0 + (-uy - uz) * (4.5 * (-uy - uz) + 3.0) - u2);
        LOCAL(srcGrid, ET) = DFL3 * rho * (1.0 + (+ux + uz) * (4.5 * (+ux + uz) + 3.0) - u2);
        LOCAL(srcGrid, EB) = DFL3 * rho * (1.0 + (+ux - uz) * (4.5 * (+ux - uz) + 3.0) - u2);
        LOCAL(srcGrid, WT) = DFL3 * rho * (1.0 + (-ux + uz) * (4.5 * (-ux + uz) + 3.0) - u2);
        LOCAL(srcGrid, WB) = DFL3 * rho * (1.0 + (-ux - uz) * (4.5 * (-ux - uz) + 3.0) - u2);
    SWEEP_END
}

void LBM_showGridStatistics(LBM_Grid grid) {
    int nObstacleCells = 0,
            nAccelCells = 0,
            nFluidCells = 0;
    float ux, uy, uz;
    float minU2 = 1e+30, maxU2 = -1e+30, u2;
    float minRho = 1e+30, maxRho = -1e+30, rho;
    float mass = 0;
    SWEEP_VAR
    SWEEP_START(0, 0, 0, 0, 0, SIZE_Z)
        rho = +LOCAL(grid, C) + LOCAL(grid, N)
              + LOCAL(grid, S) + LOCAL(grid, E)
              + LOCAL(grid, W) + LOCAL(grid, T)
              + LOCAL(grid, B) + LOCAL(grid, NE)
              + LOCAL(grid, NW) + LOCAL(grid, SE)
              + LOCAL(grid, SW) + LOCAL(grid, NT)
              + LOCAL(grid, NB) + LOCAL(grid, ST)
              + LOCAL(grid, SB) + LOCAL(grid, ET)
              + LOCAL(grid, EB) + LOCAL(grid, WT)
              + LOCAL(grid, WB);
        if (rho < minRho) minRho = rho;
        if (rho > maxRho) maxRho = rho;
        mass += rho;
        if (TEST_FLAG_SWEEP(grid, OBSTACLE)) {
            nObstacleCells++;
        } else {
            if (TEST_FLAG_SWEEP(grid, ACCEL))
                nAccelCells++;
            else
                nFluidCells++;
            ux = +LOCAL(grid, E) - LOCAL(grid, W) + LOCAL(grid, NE) - LOCAL(grid, NW) + LOCAL(grid, SE) -
                 LOCAL(grid, SW) + LOCAL(grid, ET) + LOCAL(grid, EB) - LOCAL(grid, WT) - LOCAL(grid, WB);
            uy = +LOCAL(grid, N) - LOCAL(grid, S) + LOCAL(grid, NE) + LOCAL(grid, NW) - LOCAL(grid, SE) -
                 LOCAL(grid, SW)
                 + LOCAL(grid, NT) + LOCAL(grid, NB) - LOCAL(grid, ST) - LOCAL(grid, SB);
            uz = +LOCAL(grid, T) - LOCAL(grid, B) + LOCAL(grid, NT) - LOCAL(grid, NB) + LOCAL(grid, ST) -
                 LOCAL(grid, SB) + LOCAL(grid, ET) - LOCAL(grid, EB) + LOCAL(grid, WT) - LOCAL(grid, WB);
            u2 = (ux * ux + uy * uy + uz * uz) / (rho * rho);
            if (u2 < minU2) minU2 = u2;
            if (u2 > maxU2) maxU2 = u2;
        }
    SWEEP_END
    printf("LBM_showGridStatistics:\n"
           "\tnObstacleCells: %7i nAccelCells: %7i nFluidCells: %7i\n"
           "\tminRho: %8.4f maxRho: %8.4f mass: %e\n"
           "\tminU: %e maxU: %e\n\n",
           nObstacleCells, nAccelCells, nFluidCells,
           minRho, maxRho, mass,
           sqrt(minU2), sqrt(maxU2));
}

static void storeValue(FILE *file, OUTPUT_PRECISION *v) {
    const int litteBigEndianTest = 1;
    if ((*((unsigned char *) &litteBigEndianTest)) == 0) {         /* big endian */
        const char *vPtr = (char *) v;
        char buffer[sizeof(OUTPUT_PRECISION)];
        int i;
        for (i = 0; i < sizeof(OUTPUT_PRECISION); i++)
            buffer[i] = vPtr[sizeof(OUTPUT_PRECISION) - i - 1];
        fwrite(buffer, sizeof(OUTPUT_PRECISION), 1, file);
    } else {                                                     /* little endian */
        fwrite(v, sizeof(OUTPUT_PRECISION), 1, file);
    }
}

static void loadValue(FILE *file, OUTPUT_PRECISION *v) {
    const int litteBigEndianTest = 1;
    if ((*((unsigned char *) &litteBigEndianTest)) == 0) {         /* big endian */
        char *vPtr = (char *) v;
        char buffer[sizeof(OUTPUT_PRECISION)];
        int i;
        fread(buffer, sizeof(OUTPUT_PRECISION), 1, file);
        for (i = 0; i < sizeof(OUTPUT_PRECISION); i++)
            vPtr[i] = buffer[sizeof(OUTPUT_PRECISION) - i - 1];
    } else {                                                     /* little endian */
        fread(v, sizeof(OUTPUT_PRECISION), 1, file);
    }
}

void LBM_storeVelocityField(LBM_Grid grid, const char *filename, const int binary) {
    int x, y, z;
    OUTPUT_PRECISION rho, ux, uy, uz;
    FILE *file = fopen(filename, (binary ? "wb" : "w"));
    for (z = 0; z < SIZE_Z; z++) {
        for (y = 0; y < SIZE_Y; y++) {
            for (x = 0; x < SIZE_X; x++) {
                rho = +GRID_ENTRY(grid, x, y, z, C) + GRID_ENTRY(grid, x, y, z, N) + GRID_ENTRY(grid, x, y, z, S) +
                      GRID_ENTRY(grid, x, y, z, E) + GRID_ENTRY(grid, x, y, z, W) + GRID_ENTRY(grid, x, y, z, T) +
                      GRID_ENTRY(grid, x, y, z, B) + GRID_ENTRY(grid, x, y, z, NE) + GRID_ENTRY(grid, x, y, z, NW) +
                      GRID_ENTRY(grid, x, y, z, SE) + GRID_ENTRY(grid, x, y, z, SW) + GRID_ENTRY(grid, x, y, z, NT) +
                      GRID_ENTRY(grid, x, y, z, NB) + GRID_ENTRY(grid, x, y, z, ST) + GRID_ENTRY(grid, x, y, z, SB) +
                      GRID_ENTRY(grid, x, y, z, ET)
                      + GRID_ENTRY(grid, x, y, z, EB) + GRID_ENTRY(grid, x, y, z, WT) + GRID_ENTRY(grid, x, y, z, WB);
                ux = +GRID_ENTRY(grid, x, y, z, E) - GRID_ENTRY(grid, x, y, z, W)
                     + GRID_ENTRY(grid, x, y, z, NE) - GRID_ENTRY(grid, x, y, z, NW)
                     + GRID_ENTRY(grid, x, y, z, SE) - GRID_ENTRY(grid, x, y, z, SW)
                     + GRID_ENTRY(grid, x, y, z, ET) + GRID_ENTRY(grid, x, y, z, EB)
                     - GRID_ENTRY(grid, x, y, z, WT) - GRID_ENTRY(grid, x, y, z, WB);
                uy = +GRID_ENTRY(grid, x, y, z, N) - GRID_ENTRY(grid, x, y, z, S)
                     + GRID_ENTRY(grid, x, y, z, NE) + GRID_ENTRY(grid, x, y, z, NW)
                     - GRID_ENTRY(grid, x, y, z, SE) - GRID_ENTRY(grid, x, y, z, SW)
                     + GRID_ENTRY(grid, x, y, z, NT) + GRID_ENTRY(grid, x, y, z, NB)
                     - GRID_ENTRY(grid, x, y, z, ST) - GRID_ENTRY(grid, x, y, z, SB);
                uz = +GRID_ENTRY(grid, x, y, z, T) - GRID_ENTRY(grid, x, y, z, B)
                     + GRID_ENTRY(grid, x, y, z, NT) - GRID_ENTRY(grid, x, y, z, NB)
                     + GRID_ENTRY(grid, x, y, z, ST) - GRID_ENTRY(grid, x, y, z, SB)
                     + GRID_ENTRY(grid, x, y, z, ET) - GRID_ENTRY(grid, x, y, z, EB)
                     + GRID_ENTRY(grid, x, y, z, WT) - GRID_ENTRY(grid, x, y, z, WB);
                ux /= rho;
                uy /= rho;
                uz /= rho;

                if (binary) {
                    storeValue(file, &ux);
                    storeValue(file, &uy);
                    storeValue(file, &uz);
                } else
                    fprintf(file, "%e %e %e\n", ux, uy, uz);
            }
        }
    }
    fclose(file);
}

void LBM_compareVelocityField(LBM_Grid grid, const char *filename, const int binary) {
    int x, y, z;
    float rho, ux, uy, uz;
    OUTPUT_PRECISION fileUx, fileUy, fileUz, dUx, dUy, dUz, diff2, maxDiff2 = -1e+30;
    FILE *file = fopen(filename, (binary ? "rb" : "r"));
    for (z = 0; z < SIZE_Z; z++) {
        for (y = 0; y < SIZE_Y; y++) {
            for (x = 0; x < SIZE_X; x++) {
                rho = +GRID_ENTRY(grid, x, y, z, C) + GRID_ENTRY(grid, x, y, z, N)
                      + GRID_ENTRY(grid, x, y, z, S) + GRID_ENTRY(grid, x, y, z, E)
                      + GRID_ENTRY(grid, x, y, z, W) + GRID_ENTRY(grid, x, y, z, T)
                      + GRID_ENTRY(grid, x, y, z, B) + GRID_ENTRY(grid, x, y, z, NE)
                      + GRID_ENTRY(grid, x, y, z, NW) + GRID_ENTRY(grid, x, y, z, SE)
                      + GRID_ENTRY(grid, x, y, z, SW) + GRID_ENTRY(grid, x, y, z, NT)
                      + GRID_ENTRY(grid, x, y, z, NB) + GRID_ENTRY(grid, x, y, z, ST)
                      + GRID_ENTRY(grid, x, y, z, SB) + GRID_ENTRY(grid, x, y, z, ET)
                      + GRID_ENTRY(grid, x, y, z, EB) + GRID_ENTRY(grid, x, y, z, WT)
                      + GRID_ENTRY(grid, x, y, z, WB);
                ux = +GRID_ENTRY(grid, x, y, z, E) - GRID_ENTRY(grid, x, y, z, W)
                     + GRID_ENTRY(grid, x, y, z, NE) - GRID_ENTRY(grid, x, y, z, NW)
                     + GRID_ENTRY(grid, x, y, z, SE) - GRID_ENTRY(grid, x, y, z, SW)
                     + GRID_ENTRY(grid, x, y, z, ET) + GRID_ENTRY(grid, x, y, z, EB)
                     - GRID_ENTRY(grid, x, y, z, WT) - GRID_ENTRY(grid, x, y, z, WB);
                uy = +GRID_ENTRY(grid, x, y, z, N) - GRID_ENTRY(grid, x, y, z, S)
                     + GRID_ENTRY(grid, x, y, z, NE) + GRID_ENTRY(grid, x, y, z, NW)
                     - GRID_ENTRY(grid, x, y, z, SE) - GRID_ENTRY(grid, x, y, z, SW)
                     + GRID_ENTRY(grid, x, y, z, NT) + GRID_ENTRY(grid, x, y, z, NB)
                     - GRID_ENTRY(grid, x, y, z, ST) - GRID_ENTRY(grid, x, y, z, SB);
                uz = +GRID_ENTRY(grid, x, y, z, T) - GRID_ENTRY(grid, x, y, z, B)
                     + GRID_ENTRY(grid, x, y, z, NT) - GRID_ENTRY(grid, x, y, z, NB)
                     + GRID_ENTRY(grid, x, y, z, ST) - GRID_ENTRY(grid, x, y, z, SB)
                     + GRID_ENTRY(grid, x, y, z, ET) - GRID_ENTRY(grid, x, y, z, EB)
                     + GRID_ENTRY(grid, x, y, z, WT) - GRID_ENTRY(grid, x, y, z, WB);
                ux /= rho;
                uy /= rho;
                uz /= rho;
                if (binary) {
                    loadValue(file, &fileUx);
                    loadValue(file, &fileUy);
                    loadValue(file, &fileUz);
                } else {
                    if (sizeof(OUTPUT_PRECISION) == sizeof(double)) {
                        fscanf(file, "%lf %lf %lf\n", &fileUx, &fileUy, &fileUz);
                    } else {
                        fscanf(file, "%f %f %f\n", &fileUx, &fileUy, &fileUz);
                    }
                }
                dUx = ux - fileUx;
                dUy = uy - fileUy;
                dUz = uz - fileUz;
                diff2 = dUx * dUx + dUy * dUy + dUz * dUz;
                if (diff2 > maxDiff2) maxDiff2 = diff2;
            }
        }
    }

    printf("LBM_compareVelocityField: maxDiff = %e  ==>  %s\n\n",
           sqrt(maxDiff2),
           sqrt(maxDiff2) > 1e-5 ? "##### ERROR #####" : "OK");
    fclose(file);
}

static LBM_GridPtr srcGrid, dstGrid;

void MAIN_printInfo(const MAIN_Param *param) {
    const char actionString[3][32] = {"nothing", "compare", "store"};
    const char simTypeString[3][32] = {"lid-driven cavity", "channel flow"};
    printf("MAIN_printInfo:\n"
           "\tgrid size      : %i x %i x %i = %.2f * 10^6 Cells\n"
           "\tnTimeSteps     : %i\n"
           "\tresult file    : %s\n"
           "\taction         : %s\n"
           "\tsimulation type: %s\n"
           "\tobstacle file  : %s\n\n",
           SIZE_X, SIZE_Y, SIZE_Z, 1e-6 * SIZE_X * SIZE_Y * SIZE_Z,
           param->nTimeSteps, param->resultFilename,
           actionString[param->action], simTypeString[param->simType],
           (param->obstacleFilename == NULL) ? "<none>" :
           param->obstacleFilename);
}

void CPU_MAIN_initialize(const MAIN_Param *param) {
    LBM_allocateGrid((float **) &srcGrid);
    LBM_allocateGrid((float **) &dstGrid);

    LBM_initializeGrid(*srcGrid);
    LBM_initializeGrid(*dstGrid);

    if (param->obstacleFilename != NULL) {
        LBM_loadObstacleFile(*srcGrid, param->obstacleFilename);
        LBM_loadObstacleFile(*dstGrid, param->obstacleFilename);
    }

    if (param->simType == CHANNEL) {
        LBM_initializeSpecialCellsForChannel(*srcGrid);
        LBM_initializeSpecialCellsForChannel(*dstGrid);
    } else {
        LBM_initializeSpecialCellsForLDC(*srcGrid);
        LBM_initializeSpecialCellsForLDC(*dstGrid);
    }

    LBM_showGridStatistics(*srcGrid);
}

void MAIN_parseParam(const char *input, MAIN_Param *param) {
    struct stat fileStat{};
    param->nTimeSteps = N_TIME_STEPS;
    param->obstacleFilename = const_cast<char *>(input);
    if (input != nullptr) {
        if (stat(param->obstacleFilename, &fileStat) != 0) {
            printf("MAIN_parseCommandLine: cannot stat obstacle file '%s'\n",
                   param->obstacleFilename);
            exit(1);
        }
        if (fileStat.st_size != SIZE_X * SIZE_Y * SIZE_Z + (SIZE_Y + 1) * SIZE_Z) {
            printf("MAIN_parseCommandLine:\n"
                   "\tsize of file '%s' is %i bytes\n"
                   "\texpected size is %i bytes\n",
                   param->obstacleFilename, (int) fileStat.st_size,
                   SIZE_X * SIZE_Y * SIZE_Z + (SIZE_Y + 1) * SIZE_Z);
            exit(1);
        }
    } else param->obstacleFilename = nullptr;
    param->resultFilename = "lbm.out";
    param->action = STORE;
    param->simType = LDC;
}

//=============================== GPU version =============================================
static const char *PROGRAM_SOURCE[] = {
        "#define SIZE_X (120)\n",
        "#define SIZE_Y (120)\n",
        "#define SIZE_Z (150)\n",
        "#define PADDING_X (8)\n",
        "#define PADDING_Y (0)\n",
        "#define PADDING_Z (4)\n",
        "#define PADDED_X (SIZE_X+PADDING_X)\n",
        "#define PADDED_Y (SIZE_Y+PADDING_Y)\n",
        "#define PADDED_Z (SIZE_Z+PADDING_Z)\n",
        "#define TOTAL_CELLS (SIZE_X*SIZE_Y*SIZE_Z)\n",
        "#define TOTAL_PADDED_CELLS (PADDED_X*PADDED_Y*PADDED_Z)\n",
        "#define CALC_INDEX(x,y,z,e) ( e + N_CELL_ENTRIES* ((x)+(y)*PADDED_X+(z)*PADDED_X*PADDED_Y) )\n",
        "#define MARGIN (CALC_INDEX(0, 0, 2, 0) - CALC_INDEX(0,0,0,0))\n",
        "#if 1\n",
        "#define GATHER\n",
        "#else\n",
        "#define SCATTER\n",
        "#endif\n",
        "#define BLOCK_SIZE SIZE_X\n",
        "typedef enum {C = 0,\n",
        "              N, S, E, W, T, B,\n",
        "              NE, NW, SE, SW,\n",
        "              NT, NB, ST, SB,\n",
        "              ET, EB, WT, WB,\n",
        "              FLAGS, N_CELL_ENTRIES} CELL_ENTRIES;\n",
        "\n",
        "#define N_DISTR_FUNCS FLAGS\n",
        "typedef enum {OBSTACLE    = 1 << 0,\n",
        "              ACCEL       = 1 << 1,\n",
        "              IN_OUT_FLOW = 1 << 2} CELL_FLAGS;\n",
        "#define OMEGA (1.95f)\n",
        "#define OUTPUT_PRECISION float\n",
        "#define BOOL int\n",
        "#define TRUE (-1)\n",
        "#define FALSE (0)\n",
        "#define DFL1 (1.0f/ 3.0f)\n",
        "#define DFL2 (1.0f/18.0f)\n",
        "#define DFL3 (1.0f/36.0f)\n",
        "typedef float* LBM_Grid;\n",
        "typedef LBM_Grid* LBM_GridPtr;\n",
        "#define SWEEP_X  __temp_x__\n",
        "#define SWEEP_Y  __temp_y__\n",
        "#define SWEEP_Z  __temp_z__\n",
        "#define SWEEP_VAR int __temp_x__, __temp_y__, __temp_z__;\n",
        "#define SWEEP_START(x1,y1,z1,x2,y2,z2) for( __temp_z__ = z1; __temp_z__ < z2; __temp_z__++) { for( __temp_y__ = 0; __temp_y__ < SIZE_Y; __temp_y__++) { for(__temp_x__ = 0; __temp_x__ < SIZE_X; __temp_x__++) {\n",
        "#define SWEEP_END }}}\n",
        "#define GRID_ENTRY(g,x,y,z,e)          ((g)[CALC_INDEX( x,  y,  z, e)])\n",
        "#define GRID_ENTRY_SWEEP(g,dx,dy,dz,e) ((g)[CALC_INDEX((dx)+SWEEP_X, (dy)+SWEEP_Y, (dz)+SWEEP_Z, e)])\n",
        "#define LOCAL(g,e)       (GRID_ENTRY_SWEEP( g,  0,  0,  0, e ))\n",
        "#define NEIGHBOR_C(g,e)  (GRID_ENTRY_SWEEP( g,  0,  0,  0, e ))\n",
        "#define NEIGHBOR_N(g,e)  (GRID_ENTRY_SWEEP( g,  0, +1,  0, e ))\n",
        "#define NEIGHBOR_S(g,e)  (GRID_ENTRY_SWEEP( g,  0, -1,  0, e ))\n",
        "#define NEIGHBOR_E(g,e)  (GRID_ENTRY_SWEEP( g, +1,  0,  0, e ))\n",
        "#define NEIGHBOR_W(g,e)  (GRID_ENTRY_SWEEP( g, -1,  0,  0, e ))\n",
        "#define NEIGHBOR_T(g,e)  (GRID_ENTRY_SWEEP( g,  0,  0, +1, e ))\n",
        "#define NEIGHBOR_B(g,e)  (GRID_ENTRY_SWEEP( g,  0,  0, -1, e ))\n",
        "#define NEIGHBOR_NE(g,e) (GRID_ENTRY_SWEEP( g, +1, +1,  0, e ))\n",
        "#define NEIGHBOR_NW(g,e) (GRID_ENTRY_SWEEP( g, -1, +1,  0, e ))\n",
        "#define NEIGHBOR_SE(g,e) (GRID_ENTRY_SWEEP( g, +1, -1,  0, e ))\n",
        "#define NEIGHBOR_SW(g,e) (GRID_ENTRY_SWEEP( g, -1, -1,  0, e ))\n",
        "#define NEIGHBOR_NT(g,e) (GRID_ENTRY_SWEEP( g,  0, +1, +1, e ))\n",
        "#define NEIGHBOR_NB(g,e) (GRID_ENTRY_SWEEP( g,  0, +1, -1, e ))\n",
        "#define NEIGHBOR_ST(g,e) (GRID_ENTRY_SWEEP( g,  0, -1, +1, e ))\n",
        "#define NEIGHBOR_SB(g,e) (GRID_ENTRY_SWEEP( g,  0, -1, -1, e ))\n",
        "#define NEIGHBOR_ET(g,e) (GRID_ENTRY_SWEEP( g, +1,  0, +1, e ))\n",
        "#define NEIGHBOR_EB(g,e) (GRID_ENTRY_SWEEP( g, +1,  0, -1, e ))\n",
        "#define NEIGHBOR_WT(g,e) (GRID_ENTRY_SWEEP( g, -1,  0, +1, e ))\n",
        "#define NEIGHBOR_WB(g,e) (GRID_ENTRY_SWEEP( g, -1,  0, -1, e ))\n",
        "\n",
        "#ifdef SCATTER\n",
        "\n",
        "#define SRC_C(g)  (LOCAL( g, C  ))\n",
        "#define SRC_N(g)  (LOCAL( g, N  ))\n",
        "#define SRC_S(g)  (LOCAL( g, S  ))\n",
        "#define SRC_E(g)  (LOCAL( g, E  ))\n",
        "#define SRC_W(g)  (LOCAL( g, W  ))\n",
        "#define SRC_T(g)  (LOCAL( g, T  ))\n",
        "#define SRC_B(g)  (LOCAL( g, B  ))\n",
        "#define SRC_NE(g) (LOCAL( g, NE ))\n",
        "#define SRC_NW(g) (LOCAL( g, NW ))\n",
        "#define SRC_SE(g) (LOCAL( g, SE ))\n",
        "#define SRC_SW(g) (LOCAL( g, SW ))\n",
        "#define SRC_NT(g) (LOCAL( g, NT ))\n",
        "#define SRC_NB(g) (LOCAL( g, NB ))\n",
        "#define SRC_ST(g) (LOCAL( g, ST ))\n",
        "#define SRC_SB(g) (LOCAL( g, SB ))\n",
        "#define SRC_ET(g) (LOCAL( g, ET ))\n",
        "#define SRC_EB(g) (LOCAL( g, EB ))\n",
        "#define SRC_WT(g) (LOCAL( g, WT ))\n",
        "#define SRC_WB(g) (LOCAL( g, WB ))\n",
        "\n",
        "#define DST_C(g)  (NEIGHBOR_C ( g, C  ))\n",
        "#define DST_N(g)  (NEIGHBOR_N ( g, N  ))\n",
        "#define DST_S(g)  (NEIGHBOR_S ( g, S  ))\n",
        "#define DST_E(g)  (NEIGHBOR_E ( g, E  ))\n",
        "#define DST_W(g)  (NEIGHBOR_W ( g, W  ))\n",
        "#define DST_T(g)  (NEIGHBOR_T ( g, T  ))\n",
        "#define DST_B(g)  (NEIGHBOR_B ( g, B  ))\n",
        "#define DST_NE(g) (NEIGHBOR_NE( g, NE ))\n",
        "#define DST_NW(g) (NEIGHBOR_NW( g, NW ))\n",
        "#define DST_SE(g) (NEIGHBOR_SE( g, SE ))\n",
        "#define DST_SW(g) (NEIGHBOR_SW( g, SW ))\n",
        "#define DST_NT(g) (NEIGHBOR_NT( g, NT ))\n",
        "#define DST_NB(g) (NEIGHBOR_NB( g, NB ))\n",
        "#define DST_ST(g) (NEIGHBOR_ST( g, ST ))\n",
        "#define DST_SB(g) (NEIGHBOR_SB( g, SB ))\n",
        "#define DST_ET(g) (NEIGHBOR_ET( g, ET ))\n",
        "#define DST_EB(g) (NEIGHBOR_EB( g, EB ))\n",
        "#define DST_WT(g) (NEIGHBOR_WT( g, WT ))\n",
        "#define DST_WB(g) (NEIGHBOR_WB( g, WB ))\n",
        "\n",
        "#else /* GATHER */\n",
        "\n",
        "#define SRC_C(g)  (NEIGHBOR_C ( g, C  ))\n",
        "#define SRC_N(g)  (NEIGHBOR_S ( g, N  ))\n",
        "#define SRC_S(g)  (NEIGHBOR_N ( g, S  ))\n",
        "#define SRC_E(g)  (NEIGHBOR_W ( g, E  ))\n",
        "#define SRC_W(g)  (NEIGHBOR_E ( g, W  ))\n",
        "#define SRC_T(g)  (NEIGHBOR_B ( g, T  ))\n",
        "#define SRC_B(g)  (NEIGHBOR_T ( g, B  ))\n",
        "#define SRC_NE(g) (NEIGHBOR_SW( g, NE ))\n",
        "#define SRC_NW(g) (NEIGHBOR_SE( g, NW ))\n",
        "#define SRC_SE(g) (NEIGHBOR_NW( g, SE ))\n",
        "#define SRC_SW(g) (NEIGHBOR_NE( g, SW ))\n",
        "#define SRC_NT(g) (NEIGHBOR_SB( g, NT ))\n",
        "#define SRC_NB(g) (NEIGHBOR_ST( g, NB ))\n",
        "#define SRC_ST(g) (NEIGHBOR_NB( g, ST ))\n",
        "#define SRC_SB(g) (NEIGHBOR_NT( g, SB ))\n",
        "#define SRC_ET(g) (NEIGHBOR_WB( g, ET ))\n",
        "#define SRC_EB(g) (NEIGHBOR_WT( g, EB ))\n",
        "#define SRC_WT(g) (NEIGHBOR_EB( g, WT ))\n",
        "#define SRC_WB(g) (NEIGHBOR_ET( g, WB ))\n",
        "\n",
        "#define DST_C(g)  (LOCAL( g, C  ))\n",
        "#define DST_N(g)  (LOCAL( g, N  ))\n",
        "#define DST_S(g)  (LOCAL( g, S  ))\n",
        "#define DST_E(g)  (LOCAL( g, E  ))\n",
        "#define DST_W(g)  (LOCAL( g, W  ))\n",
        "#define DST_T(g)  (LOCAL( g, T  ))\n",
        "#define DST_B(g)  (LOCAL( g, B  ))\n",
        "#define DST_NE(g) (LOCAL( g, NE ))\n",
        "#define DST_NW(g) (LOCAL( g, NW ))\n",
        "#define DST_SE(g) (LOCAL( g, SE ))\n",
        "#define DST_SW(g) (LOCAL( g, SW ))\n",
        "#define DST_NT(g) (LOCAL( g, NT ))\n",
        "#define DST_NB(g) (LOCAL( g, NB ))\n",
        "#define DST_ST(g) (LOCAL( g, ST ))\n",
        "#define DST_SB(g) (LOCAL( g, SB ))\n",
        "#define DST_ET(g) (LOCAL( g, ET ))\n",
        "#define DST_EB(g) (LOCAL( g, EB ))\n",
        "#define DST_WT(g) (LOCAL( g, WT ))\n",
        "#define DST_WB(g) (LOCAL( g, WB ))\n",
        "\n",
        "#endif /* GATHER */\n",
        "#define MAGIC_CAST(v) ((unsigned int*) ((void*) (&(v))))\n",
        "#define FLAG_VAR(v) unsigned int* _aux_ = MAGIC_CAST(v)\n",
        "#define TEST_FLAG_SWEEP(g,f)     ((*MAGIC_CAST(LOCAL(g, FLAGS))) & (f))\n",
        "#define SET_FLAG_SWEEP(g,f)      {FLAG_VAR(LOCAL(g, FLAGS)); (*_aux_) |=  (f);}\n",
        "#define CLEAR_FLAG_SWEEP(g,f)    {FLAG_VAR(LOCAL(g, FLAGS)); (*_aux_) &= ~(f);}\n",
        "#define CLEAR_ALL_FLAGS_SWEEP(g) {FLAG_VAR(LOCAL(g, FLAGS)); (*_aux_)  =    0;}\n",
        "#define TEST_FLAG(g,x,y,z,f)     ((*MAGIC_CAST(GRID_ENTRY(g, x, y, z, FLAGS))) & (f))\n",
        "#define SET_FLAG(g,x,y,z,f)      {FLAG_VAR(GRID_ENTRY(g, x, y, z, FLAGS)); (*_aux_) |=  (f);}\n",
        "#define CLEAR_FLAG(g,x,y,z,f)    {FLAG_VAR(GRID_ENTRY(g, x, y, z, FLAGS)); (*_aux_) &= ~(f);}\n",
        "#define CLEAR_ALL_FLAGS(g,x,y,z) {FLAG_VAR(GRID_ENTRY(g, x, y, z, FLAGS)); (*_aux_)  =    0;}\n",
        "\n",
        "__kernel void performStreamCollide_kernel( __global float* srcGrid, __global float* dstGrid ){\n",
        "	srcGrid += MARGIN;\n",
        "	dstGrid += MARGIN;\n",
        "    SWEEP_VAR\n",
        "	SWEEP_X = get_local_id(0);\n",
        "	SWEEP_Y = get_group_id(0);\n",
        "	SWEEP_Z = get_group_id(1);\n",
        "	float temp_swp, tempC, tempN, tempS, tempE, tempW, tempT, tempB;\n",
        "	float tempNE, tempNW, tempSE, tempSW, tempNT, tempNB, tempST ;\n",
        "	float tempSB, tempET, tempEB, tempWT, tempWB ;\n",
        "	tempC = SRC_C(srcGrid);\n",
        "	tempN = SRC_N(srcGrid);\n",
        "	tempS = SRC_S(srcGrid);\n",
        "	tempE = SRC_E(srcGrid);\n",
        "	tempW = SRC_W(srcGrid);\n",
        "	tempT = SRC_T(srcGrid);\n",
        "	tempB = SRC_B(srcGrid);\n",
        "	tempNE = SRC_NE(srcGrid);\n",
        "	tempNW = SRC_NW(srcGrid);\n",
        "	tempSE = SRC_SE(srcGrid);\n",
        "	tempSW = SRC_SW(srcGrid);\n",
        "	tempNT = SRC_NT(srcGrid);\n",
        "	tempNB = SRC_NB(srcGrid);\n",
        "	tempST = SRC_ST(srcGrid);\n",
        "	tempSB = SRC_SB(srcGrid);\n",
        "	tempET = SRC_ET(srcGrid);\n",
        "	tempEB = SRC_EB(srcGrid);\n",
        "	tempWT = SRC_WT(srcGrid);\n",
        "	tempWB = SRC_WB(srcGrid);\n",
        "	if(as_uint(LOCAL(srcGrid,FLAGS)) & (OBSTACLE)) {\n",
        "		temp_swp = tempN ; tempN = tempS ; tempS = temp_swp ;\n",
        "		temp_swp = tempE ; tempE = tempW ; tempW = temp_swp;\n",
        "		temp_swp = tempT ; tempT = tempB ; tempB = temp_swp;\n",
        "		temp_swp = tempNE; tempNE = tempSW ; tempSW = temp_swp;\n",
        "		temp_swp = tempNW; tempNW = tempSE ; tempSE = temp_swp;\n",
        "		temp_swp = tempNT ; tempNT = tempSB ; tempSB = temp_swp;\n",
        "		temp_swp = tempNB ; tempNB = tempST ; tempST = temp_swp;\n",
        "		temp_swp = tempET ; tempET= tempWB ; tempWB = temp_swp;\n",
        "		temp_swp = tempEB ; tempEB = tempWT ; tempWT = temp_swp;\n",
        "	}\n",
        "	else {\n",
        "	    float ux, uy, uz, rho, u2;\n",
        "		float temp1, temp2, temp_base;\n",
        "		rho = tempC + tempN\n",
        "			+ tempS + tempE\n",
        "			+ tempW + tempT\n",
        "			+ tempB + tempNE\n",
        "			+ tempNW + tempSE\n",
        "			+ tempSW + tempNT\n",
        "			+ tempNB + tempST\n",
        "			+ tempSB + tempET\n",
        "			+ tempEB + tempWT\n",
        "			+ tempWB;\n",
        "\n",
        "		ux = + tempE - tempW\n",
        "			+ tempNE - tempNW\n",
        "			+ tempSE - tempSW\n",
        "			+ tempET + tempEB\n",
        "			- tempWT - tempWB;\n",
        "\n",
        "		uy = + tempN - tempS\n",
        "			+ tempNE + tempNW\n",
        "			- tempSE - tempSW\n",
        "			+ tempNT + tempNB\n",
        "			- tempST - tempSB;\n",
        "\n",
        "		uz = + tempT - tempB\n",
        "			+ tempNT - tempNB\n",
        "			+ tempST - tempSB\n",
        "			+ tempET - tempEB\n",
        "			+ tempWT - tempWB;\n",
        "\n",
        "		ux /= rho;\n",
        "		uy /= rho;\n",
        "		uz /= rho;\n",
        "		if(as_uint(LOCAL(srcGrid,FLAGS)) & (ACCEL)) {\n",
        "\n",
        "			ux = 0.005f;\n",
        "			uy = 0.002f;\n",
        "			uz = 0.000f;\n",
        "		}\n",
        "		u2 = 1.5f * (ux*ux + uy*uy + uz*uz) - 1.0f;\n",
        "		temp_base = OMEGA*rho;\n",
        "		temp1 = DFL1*temp_base;\n",
        "		temp_base = OMEGA*rho;\n",
        "		temp1 = DFL1*temp_base;\n",
        "		temp2 = 1.0f-OMEGA;\n",
        "		tempC = temp2*tempC + temp1*(                                 - u2);\n",
        "	        temp1 = DFL2*temp_base;\n",
        "		tempN = temp2*tempN + temp1*(       uy*(4.5f*uy       + 3.0f) - u2);\n",
        "		tempS = temp2*tempS + temp1*(       uy*(4.5f*uy       - 3.0f) - u2);\n",
        "		tempT = temp2*tempT + temp1*(       uz*(4.5f*uz       + 3.0f) - u2);\n",
        "		tempB = temp2*tempB + temp1*(       uz*(4.5f*uz       - 3.0f) - u2);\n",
        "		tempE = temp2*tempE + temp1*(       ux*(4.5f*ux       + 3.0f) - u2);\n",
        "		tempW = temp2*tempW + temp1*(       ux*(4.5f*ux       - 3.0f) - u2);\n",
        "		temp1 = DFL3*temp_base;\n",
        "		tempNT= temp2*tempNT + temp1 *( (+uy+uz)*(4.5f*(+uy+uz) + 3.0f) - u2);\n",
        "		tempNB= temp2*tempNB + temp1 *( (+uy-uz)*(4.5f*(+uy-uz) + 3.0f) - u2);\n",
        "		tempST= temp2*tempST + temp1 *( (-uy+uz)*(4.5f*(-uy+uz) + 3.0f) - u2);\n",
        "		tempSB= temp2*tempSB + temp1 *( (-uy-uz)*(4.5f*(-uy-uz) + 3.0f) - u2);\n",
        "		tempNE = temp2*tempNE + temp1 *( (+ux+uy)*(4.5f*(+ux+uy) + 3.0f) - u2);\n",
        "		tempSE = temp2*tempSE + temp1 *((+ux-uy)*(4.5f*(+ux-uy) + 3.0f) - u2);\n",
        "		tempET = temp2*tempET + temp1 *( (+ux+uz)*(4.5f*(+ux+uz) + 3.0f) - u2);\n",
        "		tempEB = temp2*tempEB + temp1 *( (+ux-uz)*(4.5f*(+ux-uz) + 3.0f) - u2);\n",
        "		tempNW = temp2*tempNW + temp1 *( (-ux+uy)*(4.5f*(-ux+uy) + 3.0f) - u2);\n",
        "		tempSW = temp2*tempSW + temp1 *( (-ux-uy)*(4.5f*(-ux-uy) + 3.0f) - u2);\n",
        "		tempWT = temp2*tempWT + temp1 *( (-ux+uz)*(4.5f*(-ux+uz) + 3.0f) - u2);\n",
        "		tempWB = temp2*tempWB + temp1 *( (-ux-uz)*(4.5f*(-ux-uz) + 3.0f) - u2);\n",
        "	}\n",
        "	DST_C ( dstGrid ) = tempC;\n",
        "	DST_N ( dstGrid ) = tempN;\n",
        "	DST_S ( dstGrid ) = tempS;\n",
        "	DST_E ( dstGrid ) = tempE;\n",
        "	DST_W ( dstGrid ) = tempW;\n",
        "	DST_T ( dstGrid ) = tempT;\n",
        "	DST_B ( dstGrid ) = tempB;\n",
        "	DST_NE( dstGrid ) = tempNE;\n",
        "	DST_NW( dstGrid ) = tempNW;\n",
        "	DST_SE( dstGrid ) = tempSE;\n",
        "	DST_SW( dstGrid ) = tempSW;\n",
        "	DST_NT( dstGrid ) = tempNT;\n",
        "	DST_NB( dstGrid ) = tempNB;\n",
        "	DST_ST( dstGrid ) = tempST;\n",
        "	DST_SB( dstGrid ) = tempSB;\n",
        "	DST_ET( dstGrid ) = tempET;\n",
        "	DST_EB( dstGrid ) = tempEB;\n",
        "	DST_WT( dstGrid ) = tempWT;\n",
        "	DST_WB( dstGrid ) = tempWB;\n",
        "}\n"
};

//Unchangeable settings: volume simulation size for the given example
#define SIZE_X (120)
#define SIZE_Y (120)
#define SIZE_Z (150)

//Padding in each dimension
#define PADDING_X (8)
#define PADDING_Y (0)
#define PADDING_Z (4)

//Pitch in each dimension
#define PADDED_X (SIZE_X+PADDING_X)
#define PADDED_Y (SIZE_Y+PADDING_Y)
#define PADDED_Z (SIZE_Z+PADDING_Z)
#define TOTAL_PADDED_CELLS (PADDED_X*PADDED_Y*PADDED_Z)
#define CALC_INDEX(x, y, z, e) ( e + N_CELL_ENTRIES* ((x)+(y)*PADDED_X+(z)*PADDED_X*PADDED_Y) )
#define MARGIN (CALC_INDEX(0, 0, 2, 0) - CALC_INDEX(0,0,0,0))

#if 1
#define GATHER
#else
#define SCATTER
#endif

//OpenCL block size (not trivially changeable here)
#define BLOCK_SIZE SIZE_X
static const cl_uint PROGRAM_SOURCE_LEN = sizeof(PROGRAM_SOURCE) / sizeof(const char *);
cl_wrapper *cl;
cl_program program;
cl_command_queue exeCmdQueue;
cl_kernel kernel;

void GLBM_allocateGrid(ZeroCopyMem<float> *ptr) {
    const size_t size = TOTAL_PADDED_CELLS * N_CELL_ENTRIES * sizeof(float);
    *ptr = init_zero_copy_region_byte<float>(cl, size);
    memset(ptr->hostPtr, 0, size);
    printf("LBM_allocateGrid: allocated %.1f MByte\n", size / (1024.0 * 1024.0));
    ptr->hostPtr += MARGIN;
}


void initCL() {
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    exeCmdQueue = cl->createProfilingCmdQueue();
    program = cl->createProgram(PROGRAM_SOURCE, PROGRAM_SOURCE_LEN);
    kernel = cl->createKernel(program, "performStreamCollide_kernel");
}

void OpenCL_LBM_performStreamCollide(cl_mem *srcGrid, cl_mem *dstGrid) {
    cl_int clStatus;
    clStatus = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) srcGrid);
    OCL_ERRCK_RETVAL(clStatus);
    clStatus = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) dstGrid);
    OCL_ERRCK_RETVAL(clStatus);
    size_t dimBlock[3] = {SIZE_X, 1, 1};
    size_t dimGrid[3] = {SIZE_X * SIZE_Y, SIZE_Z, 1};
    cl_event event;
    clStatus = clEnqueueNDRangeKernel(exeCmdQueue, kernel, 3, nullptr, dimGrid, dimBlock, 0, NULL, &event);
    clWaitForEvents(1, &event);
    OCL_ERRCK_RETVAL(clStatus);
}

void LBM_swapGrids(ZeroCopyMem<float> *grid1, ZeroCopyMem<float> *grid2) {
    ZeroCopyMem<float> *aux = grid1;
    grid1 = grid2;
    grid2 = aux;
}

#endif //MOBILEHETEROGENOUSPROJECT_CLION_LBM_H
