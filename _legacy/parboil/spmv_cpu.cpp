//
// Created by pfxu on 7/20/19.
//
#include "heteroCompLib.h"
#include <endian.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <inttypes.h>

typedef struct _mat_entry {
    int row, col; /* i,j */
    float val;
} mat_entry;

typedef struct _row_stats {
    int index;
    int size;
    int start;
    int padding;
} row_stats;

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

typedef char MM_typecode[4];

char *mm_typecode_to_str(MM_typecode matcode);

int mm_read_banner(FILE *f, MM_typecode *matcode);

int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz);

/********************* MM_typecode query fucntions ***************************/
#define mm_is_matrix(typecode)    ((typecode)[0]=='M')
#define mm_is_sparse(typecode)    ((typecode)[1]=='C')
#define mm_is_dense(typecode)    ((typecode)[1]=='A')
#define mm_is_complex(typecode)    ((typecode)[2]=='C')
#define mm_is_real(typecode)        ((typecode)[2]=='R')
#define mm_is_pattern(typecode)    ((typecode)[2]=='P')
#define mm_is_integer(typecode) ((typecode)[2]=='I')
#define mm_is_symmetric(typecode)((typecode)[3]=='S')
#define mm_is_general(typecode)    ((typecode)[3]=='G')
#define mm_is_skew(typecode)    ((typecode)[3]=='K')
#define mm_is_hermitian(typecode)((typecode)[3]=='H')

/********************* MM_typecode modify fucntions ***************************/
#define mm_set_matrix(typecode)    ((*typecode)[0]='M')
#define mm_set_coordinate(typecode)    ((*typecode)[1]='C')
#define mm_set_array(typecode)    ((*typecode)[1]='A')
#define mm_set_dense(typecode)    mm_set_array(typecode)
#define mm_set_sparse(typecode)    mm_set_coordinate(typecode)
#define mm_set_complex(typecode)((*typecode)[2]='C')
#define mm_set_real(typecode)    ((*typecode)[2]='R')
#define mm_set_pattern(typecode)((*typecode)[2]='P')
#define mm_set_integer(typecode)((*typecode)[2]='I')
#define mm_set_symmetric(typecode)((*typecode)[3]='S')
#define mm_set_general(typecode)((*typecode)[3]='G')
#define mm_set_skew(typecode)    ((*typecode)[3]='K')
#define mm_set_hermitian(typecode)((*typecode)[3]='H')
#define mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= (*typecode)[2]=' ',(*typecode)[3]='G')
#define MM_PREMATURE_EOF        12
#define MM_NO_HEADER            14
#define MM_UNSUPPORTED_TYPE        15
#define MM_MTX_STR        "matrix"
#define MM_DENSE_STR    "array"
#define MM_SPARSE_STR    "coordinate"
#define MM_COMPLEX_STR    "complex"
#define MM_REAL_STR        "real"
#define MM_INT_STR        "integer"
#define MM_GENERAL_STR  "general"
#define MM_SYMM_STR        "symmetric"
#define MM_HERM_STR        "hermitian"
#define MM_SKEW_STR        "skew-symmetric"
#define MM_PATTERN_STR  "pattern"


int sort_rows(const void *a, const void *b) {
    return (((mat_entry *) a)->row - ((mat_entry *) b)->row);
}

/* sorts largest by size first */
int sort_stats(const void *a, const void *b) {
    return (((row_stats *) b)->size - ((row_stats *) a)->size);
}

#if __BYTE_ORDER != __LITTLE_ENDIAN
# error "File I/O is not implemented for this system: wrong endianness."
#endif

int mm_read_banner(FILE *f, MM_typecode *matcode) {
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH];
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;

    mm_clear_typecode(matcode);

    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
        return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type,
               storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    for (p = mtx; *p != '\0'; *p = tolower(*p), p++);  /* convert to lower case */
    for (p = crd; *p != '\0'; *p = tolower(*p), p++);
    for (p = data_type; *p != '\0'; *p = tolower(*p), p++);
    for (p = storage_scheme; *p != '\0'; *p = tolower(*p), p++);

    /* check for banner */
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    /* first field should be "mtx" */
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return MM_UNSUPPORTED_TYPE;
    mm_set_matrix(matcode);

    if (strcmp(crd, MM_SPARSE_STR) == 0)
        mm_set_sparse(matcode);
    else if (strcmp(crd, MM_DENSE_STR) == 0)
        mm_set_dense(matcode);
    else
        return MM_UNSUPPORTED_TYPE;

    if (strcmp(data_type, MM_REAL_STR) == 0)
        mm_set_real(matcode);
    else if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        mm_set_complex(matcode);
    else if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
    else if (strcmp(data_type, MM_INT_STR) == 0)
        mm_set_integer(matcode);
    else
        return MM_UNSUPPORTED_TYPE;

    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        mm_set_general(matcode);
    else if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        mm_set_symmetric(matcode);
    else if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
    else if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        mm_set_skew(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
    return 0;
}

int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz) {
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;
    *M = *N = *nz = 0;
    do {
        if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
            return MM_PREMATURE_EOF;
    } while (line[0] == '%');
    if (sscanf(line, "%d %d %d", M, N, nz) == 3)
        return 0;
    else
        do {
            num_items_read = fscanf(f, "%d %d %d", M, N, nz);
            if (num_items_read == EOF) return MM_PREMATURE_EOF;
        } while (num_items_read != 3);
    return 0;
}

char *mm_strdup(const char *s) {
    int len = strlen(s);
    char *s2 = (char *) malloc((len + 1) * sizeof(char));
    return strcpy(s2, s);
}

char *mm_typecode_to_str(MM_typecode matcode) {
    char buffer[MM_MAX_LINE_LENGTH];
    char *types[4];
    char *mm_strdup(const char *);
    int error = 0;
    if (mm_is_matrix(matcode))
        types[0] = MM_MTX_STR;
    else
        error = 1;
    if (mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;
    if (mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else if (mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;
    if (mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else if (mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else if (mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    sprintf(buffer, "%s %s %s %s", types[0], types[1], types[2], types[3]);
    return mm_strdup(buffer);

}


int coo_to_jds(char *mtx_filename, int pad_rows, int warp_size, int pack_size,
               int mirrored, int binary, int debug_level,
               float **data, int **data_row_ptr, int **nz_count, int **data_col_index,
               int **data_row_map, int *data_cols, int *dim, int *len, int *nz_count_len,
               int *data_ptr_len) {
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int nz;
    int i;
    float *val;
    mat_entry *entries;
    row_stats *stats;
    int rows, cols;
    if ((f = fopen(mtx_filename, "r")) == NULL)
        exit(1);
    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode)) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }
    if ((ret_code = mm_read_mtx_crd_size(f, &rows, &cols, &nz)) != 0)
        exit(1);
    *dim = rows;
    if (mirrored) {
        entries = (mat_entry *) malloc(2 * nz * sizeof(mat_entry));
    } else {
        entries = (mat_entry *) malloc(nz * sizeof(mat_entry));
    }
    int cur_i = 0;
    for (i = 0; i < nz; i++, cur_i++) {
        if (!binary) {
            fscanf(f, "%d %d %f\n", &entries[cur_i].row, &entries[cur_i].col, &entries[cur_i].val);
        } else {
            fscanf(f, "%d %d\n", &entries[cur_i].row, &entries[cur_i].col);
            entries[cur_i].val = 1.0;
        }
        entries[cur_i].row--;
        entries[cur_i].col--;
        if (mirrored) {
            if (entries[cur_i].row != entries[cur_i].col) { // not a diagonal value
                cur_i++;
                entries[cur_i].val = entries[cur_i - 1].val;
                entries[cur_i].col = entries[cur_i - 1].row;
                entries[cur_i].row = entries[cur_i - 1].col;
            }
        }
    }
    nz = cur_i;
    if (debug_level >= 1) {
        printf("Converting COO to JDS format (%dx%d)\n%d matrix entries, warp size = %d, "
               "row padding align = %d, pack size = %d\n\n", rows, cols, nz, warp_size, pad_rows, pack_size);
    }
    if (f != stdin) fclose(f);
    int irow, icol = 0, istart = 0;
    int total_size = 0;
    qsort(entries, nz, sizeof(mat_entry), sort_rows); // sort by row number
    rows = entries[nz - 1].row + 1; // last item is greatest row (zero indexed)
    if (rows % warp_size) { // pad group number to warp_size here
        rows += warp_size - rows % warp_size;
    }
    stats = (row_stats *) calloc(rows, sizeof(row_stats)); // set to 0
    *data_row_map = (int *) calloc(rows, sizeof(int));
    irow = entries[0].row; // set first row
    for (i = 0; i < nz; i++) { // loop through each sorted entry
        if (entries[i].row != irow || i == nz - 1) { // new row
            if (i == nz - 1) { icol++; }
            stats[irow].size = icol; // record # cols in previous row
            stats[irow].index = entries[i - 1].row; // row # for previous stat item
            stats[irow].start = istart; // starting location in entries array
            icol = 0; // reset row items
            irow = entries[i].row;
            istart = i;
        }
        icol++; // keep track of number of items in this row
    }
    *nz_count_len = rows / warp_size + rows % warp_size;
    *nz_count = (int *) malloc(*nz_count_len * sizeof(int)); // only one value per group
    qsort(stats, rows, sizeof(row_stats), sort_stats);
    if (debug_level >= 1) {
        printf("Padding data....%d rows, %d groups\n", rows, *nz_count_len);
    }
    int pad_to, total_padding = 0, pack_to;
    pad_rows *= pack_size; // change padding to account for packed items
    for (i = 0; i < rows; i++) {
        (*data_row_map)[i] = stats[i].index;
        if (i % warp_size == 0) { // on a group boundary with the largest number of items
            if (stats[i].size % pad_rows) {
                stats[i].padding = pad_rows - (stats[i].size % pad_rows); // find padding
            } else {
                stats[i].padding = 0; // no padding necessary, already at pad multiple
            }
            if (stats[i].size % pack_size) {
                pack_to = ceil((float) stats[i].size / pack_size);
            } else {
                pack_to = stats[i].size / pack_size;
            }
            pad_to = stats[i].size + stats[i].padding; // total size of this row, with padding
            (*nz_count)[i / warp_size] = pack_to; // number of packed items in this group
            total_size += pad_to * warp_size; // allocate size for this padded group
            if (debug_level >= 2)
                printf("Padding warp group %d to %d items, zn = %d\n", i / warp_size, pad_to, pack_to);
        } else {
            stats[i].padding = pad_to - stats[i].size;
        }
        total_padding += stats[i].padding;
    }

    if (debug_level >= 1)
        printf("Allocating data space: %d entries (%f%% padding)\n", total_size,
               (float) 100 * total_padding / total_size);
    *data = (float *) calloc(total_size, sizeof(float)); // set to 0 so padded values are set
    *data_col_index = (int *) calloc(total_size, sizeof(int)); // any unset indexes point to 0
    *data_row_ptr = (int *) calloc(rows, sizeof(int));
    *len = total_size;
    i = 0; // data index, including padding
    irow = 0; // keep track of which row we are in inside the fubini-ed array
    int idata = 0; // position within final data array
    int entry_index, j;
    int ipack; // used in internal loop for writing packed values
    mat_entry entry;
    while (1) {
        (*data_row_ptr)[irow] = idata;
        if (stats[0].size + stats[0].padding <= irow * pack_size) break;
        for (i = 0; i < rows; i++) {
            for (ipack = 0; ipack < pack_size; ipack++) {
                if (stats[i].size > irow * pack_size + ipack) {
                    entry_index = stats[i].start + irow * pack_size + ipack;
                    entry = entries[entry_index];
                    (*data)[idata] = entry.val;
                    (*data_col_index)[idata] = entry.col;
                    if (debug_level >= 2) {
                        if (i < 3) {
                            printf("[%d row%d=%.3f]", ipack + 1, i, entry.val);
                        } else {
                            printf("%d", ipack + 1);
                        }
                    }
                } else if (stats[i].size + stats[i].padding > irow * pack_size + ipack) {
                    if (debug_level >= 2) printf("0");
                    (*data_col_index)[idata] = -1;
                } else {
                    goto endwrite; // no data written this pass, so don't increment idata
                }
                idata += 1;
            }
        }
        endwrite:
        if (debug_level >= 2) {
            printf("\n");
        }
        irow += 1;
    }
    if (debug_level >= 1)
        printf("Finished converting.\nJDS format has %d columns, %d rows.\n", rows, irow);
    free(entries);
    free(stats);
    printf("nz_count_len = %d\n", *nz_count_len);
    *data_cols = rows;
    *data_ptr_len = irow + 1;
    return 0;
}

void input_vec(char *fName, float *h_vec, int dim) {
    FILE *fid = fopen(fName, "rb");
    fread(h_vec, sizeof(float), dim, fid);
    fclose(fid);
}

//void runCPU(char *input) {
    int runCPU(const char *input, double *exeTime) {
    char matPath[128];
    char vecPath[128];
    strcpy(matPath, input);
    strcpy(vecPath, input);
    strcat(matPath, "mtx.bin");
    strcat(vecPath, "vector.bin");
    //parameters declaration
    int len;
    int depth;
    int dim;
    int pad = 1;
    int nzcnt_len;
    //matrix
    float *h_data;
    int *h_indices;
    int *h_ptr;
    int *h_perm;
    int *h_nzcnt;
    //vector
    float *h_Ax_vector;
    float *h_x_vector;

    int col_count;
    coo_to_jds(matPath, 1, pad, 1, 1, 0, 1,
               &h_data, &h_ptr, &h_nzcnt, &h_indices, &h_perm,
               &col_count, &dim, &len, &nzcnt_len, &depth);

    h_Ax_vector = (float *) malloc(sizeof(float) * dim);
    h_x_vector = (float *) malloc(sizeof(float) * dim);
    input_vec(vecPath, h_x_vector, dim);

    flushed_printf("\n\tspmv Running...\n");
    double timer = getCurrentTime();
    int p, i;
    for (p = 0; p < 50; p++) {
#pragma omp parallel for
        for (i = 0; i < dim; i++) {
            int k;
            float sum = 0.0f;
            int bound = h_nzcnt[i];
            for (k = 0; k < bound; k++) {
                int j = h_ptr[k] + i;
                int in = h_indices[j];
                float d = h_data[j];
                float t = h_x_vector[in];
                sum += d * t;
            }
            //  #pragma omp critical
            h_Ax_vector[h_perm[i]] = sum;
        }
    }
    timer = getCurrentTime() - timer;
    flushed_printf("\tspmv CPU done, time: %f sec\n", timer);
    *exeTime = timer;

    free(h_data);
    free(h_indices);
    free(h_ptr);
    free(h_perm);
    free(h_nzcnt);
    free(h_Ax_vector);
    free(h_x_vector);
}

int main(int argc, char **argv) {
    char *input = "datasets/spmv/large/input/";

    std::vector<int> littleCPUList;
    std::vector<int> bigCPUList;
    for (int i = 0; i < 4; ++i) {
        littleCPUList.push_back(i);
        bigCPUList.push_back(i + 4);
    }
    char *input_small = "datasets/spmv/small/input/";
    double exetimer[8];

    flushed_printf("1 Big: \n");
    setCurThreadAffinity(7);
    omp_set_num_threads(1);
    runCPU(input, &exetimer[0]);

    flushed_printf("1 Little: \n");
    setCurThreadAffinity(3);
    omp_set_num_threads(1);
    runCPU(input, &exetimer[1]);

    flushed_printf("4 Big: \n");
    omp_set_num_threads(4);
    setCurThreadAffinity(bigCPUList);
    runCPU(input, &exetimer[2]);

    flushed_printf("4 Little: \n");
    setCurThreadAffinity(littleCPUList);
    runCPU(input, &exetimer[3]);


    flushed_printf("1 Big: \n");
    setCurThreadAffinity(7);
    omp_set_num_threads(1);
    runCPU(input_small, &exetimer[4]);

    flushed_printf("1 Little: \n");
    setCurThreadAffinity(3);
    omp_set_num_threads(1);
    runCPU(input_small, &exetimer[5]);

    flushed_printf("4 Big: \n");
    omp_set_num_threads(4);
    setCurThreadAffinity(bigCPUList);
    runCPU(input_small, &exetimer[6]);

    flushed_printf("4 Little: \n");
    setCurThreadAffinity(littleCPUList);
    runCPU(input_small, &exetimer[7]);

    printf("\n Raw data: ");
    for (int i = 0; i < 8; ++i) {
        printf("%lf, ", exetimer[i]) ;
    }
    printf("\n Ratio: ");
    for (int i = 0; i < 8; i+=2) {
        printf("%lf, ", exetimer[i + 1]/exetimer[i]) ;
    }
    printf("\n");
    return 0;
}


