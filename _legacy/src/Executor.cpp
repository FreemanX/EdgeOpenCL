//
// Created by pfxu on 2/21/19.
//

#include "Executor.h"


int Executor::executeTask(task_pack *taskPack) {
    switch (taskPack->taskType) {
        case CPU_TASK:
            
            break;

        case GPU_TASK:
            break;

        case DSP_TASK:
            // Preserve for the future
            break;
        default:
            std::cout << "Unknown task type" << std::endl;
            break;
    }
    return 0;
}
