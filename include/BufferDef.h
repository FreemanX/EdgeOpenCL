#ifndef EDCL_INCLUDE_BUFFERDEF_H_
#define EDCL_INCLUDE_BUFFERDEF_H_
#include "utils.h"
#include "Trackable.h"
/*!
 *  Define different buffer types
 *  BUF: normal opencl buffer
 *  IMG: normal opencl image
 *  Z_BUF: zero copied buffer, associated with host ptr
 *  PRIMITIVE: single primitive type value
 *  RELEASED: released buffer, don't use it
 */

// TODO: change Buffer structs to classes, determine public and private vars and funcs, add constructors and deconstructors, change buffer creation and release function.
// TODO: operator [] overload!
enum BufferType {
  BUF = 0x0a210,
//  IMG = 0x03c11,
  Z_BUF = 0x0f2b2,
  ED_BUF = 0x0e2d3,
  PRIMITIVE = 0x0ba14,
};

class InternalBuffer_ : public Trackable {
 public:
  BufferType SIGNATURE = BUF;
  size_t sizeByte{}; // external
  cl_mem mem{}; // internal
  explicit InternalBuffer_(size_t Byte) {
	SIGNATURE = BUF;
	sizeByte = Byte;
  }
  cl_mem *getCL_memPtr() { return &mem; }
};

class InternalZBuffer_ : public Trackable {
 public:
  BufferType SIGNATURE = Z_BUF;
  size_t sizeByte{}; // external
  cl_mem mem{}; // internal
  void *hostPtr = nullptr; // external
  explicit InternalZBuffer_(size_t Byte) {
	SIGNATURE = Z_BUF;
	sizeByte = Byte;
  }
  cl_mem *getCL_memPtr() { return &mem; }
};

class InternalEDCLBuffer_ : public Trackable{
 public:
  BufferType bufferType;
  void *hostPtr{};
  size_t sizeByte{};
  InternalBuffer_ *CPUBuffer{};
  InternalZBuffer_ *GPUzBuffer{};
  explicit InternalEDCLBuffer_(size_t Byte) {
	bufferType = ED_BUF;
	sizeByte = Byte;
  }
};

#define HOST_PTR(T, EDBuffer) ((T*) EDBuffer->hostPtr) // quickly get and cast hostPtr
#define FLOAT_PTR(EDBuffer) HOST_PTR(float, EDBuffer)

#endif //EDCL_INCLUDE_BUFFERDEF_H_
