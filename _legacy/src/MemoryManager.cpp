//
// Created by pfxu on 7/2/19.
//

#include "MemoryManager.hpp"
#include <math.h>
#include "utils.h"


MemoryManager::~MemoryManager() {
    for (auto &userZBuffer : this->usedZSubBufferMap)
        clReleaseMemObject(userZBuffer.second.deviceSubBuffer);

    for (auto &freeZBuffer : this->freeZSubBufferMMap)
        clReleaseMemObject(freeZBuffer.second.deviceSubBuffer);

    for (auto &zBuffer : this->zBufferList)
        clReleaseMemObject(zBuffer.deviceBuffer);

    usedZSubBufferMap.clear();
    freeZSubBufferMMap.clear();
    zBufferList.clear();
    zBufferList.shrink_to_fit();
    clReleaseEvent(clEvent);
}

MemoryManager::MemoryManager(cl_wrapper *cl) {
    this->cl = cl;
    clGetDeviceInfo(cl->device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(max_mem_alloc_size),
                    &max_mem_alloc_size, nullptr);
    clGetDeviceInfo(cl->device, CL_DEVICE_MEM_BASE_ADDR_ALIGN, sizeof(mem_base_addr_align),
                    &mem_base_addr_align, nullptr);
    DEFAULT_BUFFER_SIZE = max_mem_alloc_size / 4;
    DEFAULT_BUFFER_SIZE = 8 * 1024 * 1024;
    expendZBuffer();
    printMemStatus("MemoryManager Constructor");
}


void MemoryManager::expendZBuffer() {
    zBufferList.push_back(createZeroCopyBuffer(DEFAULT_BUFFER_SIZE));
}

ZBuffer MemoryManager::createZeroCopyBuffer(size_t size) {
    ZBuffer zBuffer{};
    zBuffer.totalSize = size;
    zBuffer.remainingSize = size;
    zBuffer.subBufferHead = 0;
    zBuffer.deviceBuffer = clCreateBuffer(
            cl->context,
            CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR,
            size, nullptr, &error);
    cl_wrapper::checkError(error);

    zBuffer.hostPtr = clEnqueueMapBuffer(
            cl->memCmdQueue,
            zBuffer.deviceBuffer, CL_TRUE, CL_MAP_WRITE, 0,
            zBuffer.totalSize, 0, nullptr, nullptr, &error);
    cl_wrapper::checkError(error);

    return zBuffer;
}

size_t MemoryManager::convert2MemAlignedSize(size_t size) {
    return (size_t) ceil(size * 1.0 / mem_base_addr_align) * mem_base_addr_align;
}

ZSubBuffer MemoryManager::createZeroCopySubBuffer(ZBuffer *zBuffer, size_t size) {
    ZSubBuffer zSubBuffer{};
    cl_buffer_region bufferRegion{.origin=zBuffer->subBufferHead, .size = size};
    zSubBuffer.size = size;
    zSubBuffer.origin = bufferRegion.origin;
    zSubBuffer.zBuffer = zBuffer;
    zSubBuffer.deviceSubBuffer = clCreateSubBuffer(
            zBuffer->deviceBuffer,
            CL_MEM_READ_WRITE,
            CL_BUFFER_CREATE_TYPE_REGION,
            &bufferRegion, &error);
    cl_wrapper::checkError(error);
    zSubBuffer.subHostPtr = clEnqueueMapBuffer(
            cl->memCmdQueue,
            zSubBuffer.deviceSubBuffer, CL_TRUE, CL_MAP_WRITE, 0,
            bufferRegion.size, 0, nullptr, &clEvent, &error);
    clWaitForEvents(1, &clEvent);
    cl_wrapper::checkError(error);

    zBuffer->subBufferHead = zBuffer->subBufferHead + size;
    zBuffer->remainingSize = zBuffer->remainingSize - size;
    return zSubBuffer;
}

void MemoryManager::insertUsedZSubBufferMap(std::pair<int, ZSubBuffer> userKeyZSubBufferPair) {
    usedZSubBufferMap.insert(userKeyZSubBufferPair);
    userKeyZSubBufferPair.second.zBuffer->numUsedSubBuffer++;
}

void MemoryManager::eraseUsedZSubBufferMap(int key) {
    usedZSubBufferMap.find(key)->second.zBuffer->numUsedSubBuffer--;
    usedZSubBufferMap.erase(key);
}

UserZBuffer MemoryManager::applyZeroCopyBuffer(size_t size) {
    printMemStatus("applyZeroCopyBuffer Called");
    size_t actualSize = convert2MemAlignedSize(size);
    UserZBuffer userZBuffer{};
    userZBuffer.key = userBufferCnt++;

// Reuse free buffer
    for (auto it = freeZSubBufferMMap.begin(); it != freeZSubBufferMMap.end();) {
        if (it->second.size > actualSize &&
            it->second.size - actualSize < SIZE_DIFF_THRESHOLD) {
            assembleUserZBuffer(&it->second, &userZBuffer, size);
            insertUsedZSubBufferMap(std::pair<int, ZSubBuffer>(userZBuffer.key, it->second));
            freeZSubBufferMMap.erase(it);
            printMemStatus("Return reused free buffer");
            printUserZBuffer(&userZBuffer);
            return userZBuffer;
        } else
            ++it;
    }

// Create sub-buffer from existing buffer if enough space
    for (auto &zBuffer : zBufferList) {
        if (zBuffer.remainingSize >= actualSize) {
            ZSubBuffer zSubBuffer = createZeroCopySubBuffer(&zBuffer, actualSize);
            assembleUserZBuffer(&zSubBuffer, &userZBuffer, size);
            insertUsedZSubBufferMap(std::pair<int, ZSubBuffer>(userZBuffer.key, zSubBuffer));
            printMemStatus("Create New Sub-buffer");
            printUserZBuffer(&userZBuffer);
            return userZBuffer;
        }
    }

// Create sub-buffer from new buffer if no enough space
    ZBuffer zBuffer = createZeroCopyBuffer(std::max(DEFAULT_BUFFER_SIZE, actualSize));
    ZSubBuffer zSubBuffer = createZeroCopySubBuffer(&zBuffer, actualSize);
    assembleUserZBuffer(&zSubBuffer, &userZBuffer, size);
    insertUsedZSubBufferMap(std::pair<int, ZSubBuffer>(userZBuffer.key, zSubBuffer));
    zBufferList.push_back(zBuffer);
    printMemStatus("Create New buffer and Sub-buffer");
    printUserZBuffer(&userZBuffer);
    return userZBuffer;
}

void MemoryManager::releaseUserZBuffer(UserZBuffer *userZBuffer) {
    ZSubBuffer zSubBuffer = usedZSubBufferMap.find(userZBuffer->key)->second;
//        ZSubBuffer zSubBuffer = usedZSubBufferMap[userZBuffer->key];
    eraseUsedZSubBufferMap(userZBuffer->key);
// if this zSubBuffer is the only sub-buffer for its zBuffer & this zBuffer isn't the only one
// then release this zSubBuffer and zBuffer it belongs to. clRelease() and remove from list&map
// else put this zSubBuffer into freeZSubBufferList
    if (zBufferList.size() > 1 && zSubBuffer.zBuffer->numUsedSubBuffer < 1) {
        freeZSubBufferMMap.erase(zSubBuffer.zBuffer);
        clReleaseMemObject(zSubBuffer.zBuffer->deviceBuffer);
//        zBufferList.erase(std::remove(zBufferList.begin(), zBufferList.end(), zSubBuffer.zBuffer),
//                          zBufferList.end());
        for (int i = 0; i < zBufferList.size(); ++i) {
            if (zBufferList[i].hostPtr == zSubBuffer.zBuffer->hostPtr
                && zBufferList[i].numUsedSubBuffer == 0) {
                zBufferList.erase(zBufferList.begin() + i);
            }
        }
    } else {
        freeZSubBufferMMap.insert(std::pair<ZBuffer *, ZSubBuffer>(zSubBuffer.zBuffer, zSubBuffer));
    }
    printMemStatus("releaseUserBuffer");
}

void MemoryManager::assembleUserZBuffer(ZSubBuffer *zSubBuffer, UserZBuffer *userZBuffer, size_t size) {
    userZBuffer->size = size;
    userZBuffer->deviceBuffer = &zSubBuffer->deviceSubBuffer;
    userZBuffer->hostPtr = zSubBuffer->subHostPtr;
}

void MemoryManager::printMemStatus(std::string msg) {
    flushed_printf(msg.c_str());
    flushed_printf("\n\tzBufferList info:\n");
    flushed_printf("\tsize: %d\n", zBufferList.size());
    for (auto zBuffer: zBufferList) {
        flushed_printf("\t\t{hostPtr: %p, ", zBuffer.hostPtr);
        flushed_printf("deviceBuffer: %p, ", &zBuffer.deviceBuffer);
        flushed_printf("totalSize: %d, ", zBuffer.totalSize);
        flushed_printf("remainingSize: %d, ", zBuffer.remainingSize);
        flushed_printf("subBufferHead: %d, ", zBuffer.subBufferHead);
        flushed_printf("numUsedSubBuffer: %d}\n", zBuffer.numUsedSubBuffer);
    }

    flushed_printf("\tusedZSubBufferMap info:\n");
    flushed_printf("\tsize: %d\n", usedZSubBufferMap.size());
    for (auto bufPair: usedZSubBufferMap) {
        flushed_printf("\t\t%d, ", bufPair.first);
        printZSubBuffer(&bufPair.second);
    }

    flushed_printf("\tfreeZSubBufferMMap info:\n");
    flushed_printf("\tsize: %d\n", freeZSubBufferMMap.size());
    for (auto bufPair: freeZSubBufferMMap) {
        flushed_printf("\t\t%p, ", bufPair.first);
        printZSubBuffer(&bufPair.second);
    }
}

void printZSubBuffer(ZSubBuffer *zSubBuffer) {
    flushed_printf("{subHostPtr: %p, ", zSubBuffer->subHostPtr);
    flushed_printf("deviceSubBuffer: %p, ", zSubBuffer->deviceSubBuffer);
    flushed_printf("*zBuffer: %p, ", zSubBuffer->zBuffer);
    flushed_printf("origin: %ld, ", zSubBuffer->origin);
    flushed_printf("size: %ld}\n", zSubBuffer->size);

}

void printUserZBuffer(UserZBuffer *userZBuffer) {
    flushed_printf("userZBuffer: {hostPtr: %p, ", userZBuffer->hostPtr);
    flushed_printf("cl_mem: %p, ", userZBuffer->deviceBuffer);
    flushed_printf("size: %ld, ", userZBuffer->size);
    flushed_printf("key: %d}\n", userZBuffer->key);
    flushed_printf("\n");
}

