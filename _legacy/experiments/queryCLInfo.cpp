#include "heteroCompLib.h"

using namespace std;

int main(int argc, char** argv)
{
    CLInfo clInfo{};
    cl_wrapper::queryCLInfo(&clInfo);

    cout << "Number of platforms: " << clInfo.num_platforms << endl;

    CLPlatform *platform;
    CLDevice *device;
    char buffer[1024];
    for(size_t i = 0; i < clInfo.num_platforms; i++)
    {
        platform = &clInfo.platforms[i];
        cl_wrapper::checkError(clGetPlatformInfo(platform->platformId, CL_PLATFORM_NAME, sizeof(buffer), &buffer, NULL));
        cout << "Platform: "  << buffer;
        cl_wrapper::checkError(clGetPlatformInfo(platform->platformId, CL_PLATFORM_VERSION, sizeof(buffer), &buffer, NULL));
        cout << "(" << buffer << "), ";
        cl_wrapper::checkError(clGetPlatformInfo(platform->platformId, CL_PLATFORM_VENDOR, sizeof(buffer), &buffer, NULL));
        cout << "Vendor: "  << buffer << endl;
        cl_wrapper::checkError(clGetPlatformInfo(platform->platformId, CL_PLATFORM_PROFILE, sizeof(buffer), &buffer, NULL));
        cout << "Platform OpenCL Profile: " << buffer << endl;
        cl_wrapper::checkError(clGetPlatformInfo(platform->platformId, CL_PLATFORM_EXTENSIONS, sizeof(buffer), &buffer, NULL));
        cout << "Supported extensions: " << buffer << endl;

        cout << "Platform " << i << " has " << platform->num_devices << " devices";
        cout << "===================" << endl;
        for(size_t j = 0; j < platform->num_devices; j++)
        {
            device  = &platform->devices[j];
            cl_wrapper::checkError(clGetDeviceInfo(device->deviceId, CL_DEVICE_NAME, sizeof(buffer), &buffer, NULL));
            cout << "Device: " << buffer;
            cl_wrapper::checkError(clGetDeviceInfo(device->deviceId, CL_DEVICE_VERSION, sizeof(buffer), &buffer, NULL));
            cout << "(" << buffer << "), ";
            cl_wrapper::checkError(clGetDeviceInfo(device->deviceId, CL_DRIVER_VERSION, sizeof(buffer), &buffer, NULL));
            cout << "software version: " << buffer << ", ";
            cl_wrapper::checkError(clGetDeviceInfo(device->deviceId, CL_DEVICE_OPENCL_C_VERSION, sizeof(buffer), &buffer, NULL));
            cout << "CL C version:" << buffer << endl;

            cout << "Number of CUs: " << device->num_CUs << endl;
            cout << "Local memory size: " << device->localMemorySizeB/1024 << "KB" << endl;
            cout << "Maximum memory allocation: " << device->maxAllocSizeB/1024 << "KB" << endl;
            cout << "WG dim size: " << device->WGDimSize << endl;
            cout << "Max WG dim sizes (";
            for (int k = 0; k < device->WGDimSize; k++)
                cout << device->maxWGSizes[k] << " ";
            cout << ")" << endl;

            cl_device_partition_property partition_property[3];
            cl_wrapper::checkError( clGetDeviceInfo(device->deviceId, CL_DEVICE_PARTITION_TYPE, sizeof(partition_property), &partition_property, nullptr) );
            cout << "Supported partition property: ";
            for (long p : partition_property) {
                if(p == CL_DEVICE_PARTITION_EQUALLY)
                    cout << "CL_DEVICE_PARTITION_EQUALLY ";
                else if(p == CL_DEVICE_PARTITION_BY_COUNTS)
                    cout << "CL_DEVICE_PARTITION_BY_COUNTS ";
                else if(p == CL_DEVICE_PARTITION_BY_AFFINITY_DOMAIN)
                    cout << "CL_DEVICE_PARTITION_BY_AFFINITY_DOMAIN ";
                else
                    cout << p << " ";
            }
            cout << endl;
        }
        cout << "===================" << endl;
    }
    
    return 0;
}