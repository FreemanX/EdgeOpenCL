//
// Created by pfxu on 7/19/19.
//

#include "benchmark.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <deque>
#include <iostream>
#include <vector>

static const char *PROGRAM_SOURCE[] = {
        "#define MAX_THREADS_PER_BLOCK 256\n",
        "#define LOCAL_MEM_SIZE 1600 //This needs to be adjusted for certain graphs with high degrees\n",
        "#define INF 2147483647//2^31-1\n",
        "#define UP_LIMIT 16677216//2^24\n",
        "#define WHITE 16677217\n",
        "#define GRAY 16677218\n",
        "#define GRAY0 16677219\n",
        "#define GRAY1 16677220\n",
        "#define BLACK 16677221\n",
        "\n",
        "struct Node {\n",
        "  int x;\n",
        "  int y;\n",
        "};\n",
        "struct Edge {\n",
        "  int x;\n",
        "  int y;\n",
        "};\n",
        "\n",
        "#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics: enable\n",
        "#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable\n",
        "#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics: enable\n",
        "__kernel void\n",
        "myKernel(__global int *q1,\n",
        "           __global int *q2,\n",
        "           __global struct Node *g_graph_nodes,\n",
        "           __global struct Edge *g_graph_edges,\n",
        "           __global int *g_color,\n",
        "           __global int *g_cost,\n",
        "           __global int *tail,\n",
        "           int no_of_nodes,\n",
        "           int gray_shade,\n",
        "           int k ,\n",
        "           __local int *local_q_tail,\n",
        "           __local int *local_q,\n",
        "           __local int *shift)\n",
        "{\n",
        "\n",
        "  if(get_local_id(0) == 0){\n",
        "    *local_q_tail = 0;//initialize the tail of w-queue\n",
        "  }\n",
        "  barrier(CLK_LOCAL_MEM_FENCE);\n",
        "\n",
        "  //first, propagate and add the new frontier elements into w-queues\n",
        "  //int tid = get_group_id(0)*MAX_THREADS_PER_BLOCK + get_local_id(0);\n",
        "  int tid = get_global_id(0);\n",
        "\n",
        "  if( tid<no_of_nodes)\n",
        "  {\n",
        "    int pid = q1[tid]; //the current frontier node, or the parent node of the new frontier nodes\n",
        "    g_color[pid] = BLACK;\n",
        "    int cur_cost = g_cost[pid];\n",
        "    //into\n",
        "    struct Node cur_node = g_graph_nodes[pid];\n",
        "\n",
        "    for(int i=cur_node.x; i<cur_node.y + cur_node.x; i++)//visit each neighbor of the\n",
        "      //current frontier node.\n",
        "    {\n",
        "      struct Edge cur_edge = g_graph_edges[i];\n",
        "      int id = cur_edge.x;\n",
        "      int cost = cur_edge.y;\n",
        "      cost += cur_cost;\n",
        "      int orig_cost = atom_min (&g_cost[id],cost);\n",
        "      if(orig_cost > cost){//the node should be visited\n",
        "        if(g_color[id] > UP_LIMIT){\n",
        "         int old_color = atom_xchg (&g_color[id],gray_shade);\n",
        "          //this guarantees that only one thread will push this node\n",
        "          //into a queue\n",
        "          if(old_color != gray_shade) {\n",
        "            //atomic operation guarantees the correctness\n",
        "            //even if multiple warps are executing simultaneously\n",
        "            int index = atom_add (local_q_tail,1);\n",
        "            local_q[index] = id;\n",
        "          }\n",
        "        }\n",
        "      }\n",
        "    }\n",
        "  }\n",
        "  barrier(CLK_LOCAL_MEM_FENCE);\n",
        "\n",
        "  if(get_local_id(0) == 0){\n",
        "    int tot_sum = *local_q_tail;\n",
        "    *shift = atom_add (tail,tot_sum);\n",
        "  }\n",
        "  barrier(CLK_LOCAL_MEM_FENCE);\n",
        "\n",
        "  //shift within a w-queue\n",
        "  int local_shift = get_local_id(0);\n",
        "\n",
        "  while(local_shift < *local_q_tail){\n",
        "    q2[*shift + local_shift] = local_q[local_shift];\n",
        "    local_shift += get_local_size(0);\n",
        "  }\n",
        "}\n",
        "\n"
};

#define MAX_THREADS_PER_BLOCK 256
#define LOCAL_MEM_SIZE 1600 //This needs to be adjusted for certain graphs with high degrees
#define INF 2147483647//2^31-1
#define UP_LIMIT 16677216//2^24
#define WHITE 16677217
#define GRAY 16677218
#define GRAY0 16677219
#define GRAY1 16677220
#define BLACK 16677221
#define INF 2147483647//2^31-1

#define WHITE 16677217
#define GRAY 16677218
#define BLACK 16677221
int no_of_nodes; //the number of nodes in the graph
int edge_list_size;//the number of edges in the graph
FILE *fp;

struct Node {
    int x;
    int y;
};

struct Edge {
    int x;
    int y;
};

////////////////////////////////////////////////////////////////////
//the cpu version of bfs for speed comparison
//the text book version ("Introduction to Algorithms")
////////////////////////////////////////////////////////////////////
void BFS_CPU(Node *h_graph_nodes, Edge *h_graph_edges, int *color, int *h_cost, int source) {
    std::deque<int> wavefront;
    wavefront.push_back(source);
    color[source] = GRAY;
    int index;
    while (!wavefront.empty()) {
        index = wavefront.front();
        wavefront.pop_front();

#pragma omp parallel for
        for (int i = h_graph_nodes[index].x;
             i < (h_graph_nodes[index].y +
                  h_graph_nodes[index].x); i++) {
            int id = h_graph_edges[i].x;
            if (color[id] == WHITE) {
                h_cost[id] = h_cost[index] + 1;

#pragma omp critical
                wavefront.push_back(id);

                color[id] = GRAY;
            }
        }
        color[index] = BLACK;
    }
}


int runCPU(const char *input, double *exeTime) {
    no_of_nodes = 0;
    edge_list_size = 0;
    flushed_printf("\tPreparing data: ");
    flushed_printf(input);
    fp = fopen(input, "r");
    if (!fp) {
        printf("\nError Reading graph file\n");
        return -1;
    }
    int source;

    fscanf(fp, "%d", &no_of_nodes);
    // allocate host memory
    Node *h_graph_nodes = (Node *) malloc(sizeof(Node) * no_of_nodes);
    int *color = (int *) malloc(sizeof(int) * no_of_nodes);
    int start, edgeno;
    // initalize the memory
    for (unsigned int i = 0; i < no_of_nodes; i++) {
        fscanf(fp, "%d %d", &start, &edgeno);
        h_graph_nodes[i].x = start;
        h_graph_nodes[i].y = edgeno;
        color[i] = WHITE;
    }
    //read the source node from the file
    fscanf(fp, "%d", &source);
    fscanf(fp, "%d", &edge_list_size);
    int id, cost;
    Edge *h_graph_edges = (Edge *) malloc(sizeof(Edge) * edge_list_size);
    for (int i = 0; i < edge_list_size; i++) {
        fscanf(fp, "%d", &id);
        fscanf(fp, "%d", &cost);
        h_graph_edges[i].x = id;
        h_graph_edges[i].y = cost;
    }
    if (fp)
        fclose(fp);

    // allocate mem for the result on host side
    int *h_cost = (int *) malloc(sizeof(int) * no_of_nodes);
    for (int i = 0; i < no_of_nodes; i++) {
        h_cost[i] = INF;
    }
    h_cost[source] = 0;

    flushed_printf("\n\tBFS Running...\n");
    double timer = getCurrentTime();
    BFS_CPU(h_graph_nodes, h_graph_edges, color, h_cost, source);
    timer = getCurrentTime() - timer;
    flushed_printf("\tBFS CPU done, time: %f sec\n", timer);

    *exeTime = timer;
    free(h_graph_nodes);
    free(color);
    free(h_cost);

    return 0;
}

static const cl_uint PROGRAM_SOURCE_LEN = sizeof(PROGRAM_SOURCE) / sizeof(const char *);
cl_wrapper *cl;
cl_program program;
cl_kernel kernel;
cl_command_queue exeCmdQueue;

void initCL() {
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);
    cl = new cl_wrapper(clInfo.platforms[0].platformId, clInfo.platforms[0].devices[0].deviceId);
    program = cl->createProgram(PROGRAM_SOURCE, PROGRAM_SOURCE_LEN);
    kernel = cl->createKernel("myKernel", program);
    exeCmdQueue = cl->createProfilingCmdQueue();
}

int runGPU(const char *input, double *exeTime) {
    initCL();
    int num_of_nodes = 0;
    int num_of_edges = 0;
    const int h_top = 1;

    fp = fopen(input, "r");
    if (!fp) {
        printf("Error Reading graph file\n");
        return 0;
    }
    int source;

    fscanf(fp, "%d", &num_of_nodes);
    // allocate host memory
    ZeroCopyMem<struct Node> graph_nodes_z = init_zero_copy_region<struct Node>(cl, num_of_nodes);
    ZeroCopyMem<int> color_z = init_zero_copy_region<int>(cl, num_of_nodes);
    int start, edgeno;
    // initialize the memory
    int i;
    for (i = 0; i < num_of_nodes; i++) {
        fscanf(fp, "%d %d", &start, &edgeno);
        graph_nodes_z.hostPtr[i].x = start;
        graph_nodes_z.hostPtr[i].y = edgeno;
        color_z.hostPtr[i] = WHITE;
    }
    //read the source node from the file
    fscanf(fp, "%d", &source);
    fscanf(fp, "%d", &num_of_edges);
    int id, cost;
    ZeroCopyMem<struct Edge> graph_edges_z = init_zero_copy_region<struct Edge>(cl, num_of_edges);
    for (i = 0; i < num_of_edges; i++) {
        fscanf(fp, "%d", &id);
        fscanf(fp, "%d", &cost);
        graph_edges_z.hostPtr[i].x = id;
        graph_edges_z.hostPtr[i].y = cost;
    }
    if (fp)
        fclose(fp);

    // allocate mem for the result on host side
    ZeroCopyMem<int> cost_z = init_zero_copy_region<int>(cl, num_of_nodes);
    for (i = 0; i < num_of_nodes; i++) { cost_z.hostPtr[i] = INF; }
    cost_z.hostPtr[source] = 0;
    ZeroCopyMem<int> q1_z = init_zero_copy_region<int>(cl, num_of_nodes);
    ZeroCopyMem<int> q2_z = init_zero_copy_region<int>(cl, num_of_nodes);
    ZeroCopyMem<int> tail_z = init_zero_copy_region<int>(cl, 1);
    tail_z.hostPtr[0] = h_top;
    cost_z.hostPtr[0] = 0;
    q1_z.hostPtr[0] = source;
    double time = getCurrentTime();
    clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &graph_nodes_z.deviceBuffer);
    clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &graph_edges_z.deviceBuffer);
    clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &color_z.deviceBuffer);
    clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &cost_z.deviceBuffer);
    clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &tail_z.deviceBuffer);
    int num_t;
    int k = 0;
    int num_of_blocks, num_of_threads_per_block;
    cl_event event;
    do {
        num_t = tail_z.hostPtr[0];
        tail_z.hostPtr[0] = 0;
        if (num_t == 0) break;
        num_of_blocks = (int) ceil(num_t / (double) MAX_THREADS_PER_BLOCK);
        num_of_threads_per_block = num_t > MAX_THREADS_PER_BLOCK ? MAX_THREADS_PER_BLOCK : num_t;

        size_t grid[1] = {static_cast<size_t>(num_of_blocks * num_of_threads_per_block)};
        size_t block[1] = {static_cast<size_t>(num_of_threads_per_block)};

        clSetKernelArg(kernel, 7, sizeof(int), (void *) &num_t);
        clSetKernelArg(kernel, 9, sizeof(int), (void *) &k);
        clSetKernelArg(kernel, 10, sizeof(int), NULL);
        clSetKernelArg(kernel, 11, LOCAL_MEM_SIZE * sizeof(int), NULL);
        clSetKernelArg(kernel, 12, sizeof(int), NULL);
        if (k % 2 == 0) {
            int gray = GRAY0;
            clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &q1_z.deviceBuffer);
            clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &q2_z.deviceBuffer);
            clSetKernelArg(kernel, 8, sizeof(int), (void *) &gray);
        } else {
            int gray = GRAY1;
            clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &q2_z.deviceBuffer);
            clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &q1_z.deviceBuffer);
            clSetKernelArg(kernel, 8, sizeof(int), (void *) &gray);
        }
        clEnqueueNDRangeKernel(exeCmdQueue, kernel, 1, 0, grid, block, 0, 0, &event);
        clWaitForEvents(1, &event);
//        printf("cl exe time: %lf\n", cl_wrapper::getExecutionTime(&event));
        k++;
//        printf("%d\n", k);
    } while (true);
    time = getCurrentTime() - time;
    flushed_printf("\tBFS GPU done, time: %f sec\n", time);

    *exeTime = time;
    return 0;
}

int main(int argc, char **argv) {
    char const *input = "datasets/bfs/SF/input/graph_input.dat";
    char const *input_small = "datasets/bfs/UT/input/graph_input.dat";
    benchmark(input, input_small, runCPU, runGPU);

    return 0;
}




















