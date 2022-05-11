//
// Created by pfxu on 4/4/19.
//

#ifndef MOBILEHETEROGENOUSPROJECT_CLION_TASKPACK_H
#define MOBILEHETEROGENOUSPROJECT_CLION_TASKPACK_H

enum taskTypes {
    CPU_TASK, GPU_TASK, DSP_TASK
};

typedef struct {
    int num_args;
    int *arg_lengths;
    void **arg_pointers;
} io_pack;

typedef struct {
    int taskType;
    int offset;
    int length;
} task_pack;

typedef struct {
    int pro_id;
    int thread_id;
    task_pack *taskPack;
    io_pack *input;
    io_pack *output;

    void (*Kernel)(io_pack *input, io_pack *output, int offset, int length);
} TaskPack;

typedef struct {
    int cpu_id;
    int thread_id;
    task_pack *taskPack;
    io_pack *input;
    io_pack *output;

    void (*cpuKernel)(io_pack *input, io_pack *output, int offset, int length);
} cpu_task_pack;


#endif //MOBILEHETEROGENOUSPROJECT_CLION_TASKPACK_H
