#!/bin/bash
ANDROID_NDK=$HOME/toolchains/android-ndk-r17c
ANDROID_TOOLCHAIN=$HOME/toolchains/ndk_r17c

ANDROID_TOOLCHAIN_BIN=$ANDROID_TOOLCHAIN/bin/
ANDROID_TOOLCHAIN_SYSROOT_USR=$ANDROID_TOOLCHAIN/sysroot/usr/

export PATH=$HOME/toolchains/platform-tools:$ANDROID_TOOLCHAIN/aarch64-linux-android/bin/:$ANDROID_TOOLCHAIN_BIN:$PATH
export LD_LIBRARY_PATH=$ANDROID_TOOLCHAIN_SYSROOT_USR/lib:$LD_LIBRARY_PATH
export HOST=aarch64-linux-android
export CC=$ANDROID_TOOLCHAIN_BIN/$HOST-clang
export CXX=$ANDROID_TOOLCHAIN_BIN/$HOST-clang++
export AR=$ANDROID_TOOLCHAIN_BIN/$HOST-ar
export LD=$ANDROID_TOOLCHAIN_BIN/$HOST-ld.gold
export RANLIB=$ANDROID_TOOLCHAIN_BIN/$HOST-ranlib
export LDFLAGS=" -pie "

reStart=0
reBuild=0
push2Phone=0
makej8=0
while [ "$1" != "" ]; do
  case $1 in
    -R | --restart)
      reStart=1
      ;;
    -r | --rebuild)
      reBuild=1
      ;;
    -p | --push)
      push2Phone=1
      ;;
    -f | --fast)
      makej8=1
      ;;
    -mp | --makePush)
      makej8=1
      push2Phone=1
      ;;
  esac
  shift
done

if [ "$reStart" = "1" ]; then
  echo "Rebuild everything"
  rm -rf build-android
fi

[ ! -d "build-android" ] && mkdir build-android
cd build-android

if [ ! "$(ls -A .)" ] || [ "$reBuild" = "1" ]; then
  cmake \
    -DANDROID_TOOLCHAIN=$ANDROID_TOOLCHAIN \
    -DCMAKE_CROSSCOMPILING=ON \
    -DCMAKE_TOOLCHAIN_FILE=../androideabi.cmake \
    -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_BUILD_TYPE:STRING=Debug \
    -DCMAKE_AR:FILEPATH=$HOST-ar \
    -DCMAKE_RANLIB:FILEPATH=$HOST-ranlib \
    -DCMAKE_CXX_FLAGS:STRING=" -ggdb -O3 -static-libstdc++ -std=c++11 -fPIE -fPIC -fuse-ld=$LD -fopenmp" \
    -DCMAKE_C_FLAGS:STRING=" -ggdb -O3 -fPIE -fPIC -fopenmp" \
    .. && ln -fs $PWD/compile_commands.json ../
fi

if [ "$makej8" = "1" ]; then
  make -j8
else
  make
fi

if [ "$push2Phone" = "1" ]; then
  adb push ./examples /data/local/tmp/EDCL
fi
