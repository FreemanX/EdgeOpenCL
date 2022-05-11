//
// Created by pfxu on 2/21/19.
// Singleton design follows:
//https://stackoverflow.com/questions/1008019/c-singleton-design-pattern
//

#ifndef PROJECT_EXECUTOR_H
#define PROJECT_EXECUTOR_H

#include <iostream>
#include "CPUExecutor.h"
#include "TaskPack.h"


class Executor {

public:
    static Executor &getInstance() {
        static Executor instance;
        return instance;
    }

    int executeTask(task_pack *taskPack);

    Executor(Executor const &) = delete;

    void operator=(Executor const &) = delete;

    ~Executor();

private:

    Executor() {};

};


#endif //PROJECT_EXECUTOR_H
