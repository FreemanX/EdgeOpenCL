#CC=arm-linux-androideabi-clang++
LIBS=-L/home/pfxu/Projects/LibOpenCL/lib -lOpenCL -llog 
INC=-I/home/pfxu/Projects/LibOpenCL/include 
CFLAGS= -g -fPIC -fPIE -static-libstdc++ -std=c++11 $(LIBS) $(INC)
LDFLAGS=-pie -fuse-ld=gold



clmempatterns: main.o
	${CXX} $(CFLAGS) $(LDFLAGS) -o $@ $^ 

main.o: main.c kernel.h
	${CXX} $(CFLAGS) $(LDFLAGS) -c $<

kernel.h: kernel.cl
	echo "const char c_kernel[]={" >kernel.h
	hexdump -ve '1/1 "0x%.2x,"' kernel.cl >>kernel.h
	echo "0x00};" >>kernel.h

clean:
	rm clmempatterns main.o kernel.h


#main: main.cpp
#	$(CXX) $(CFLAGS) $(LDFLAGS) -o cltest main.cpp 
#
#clean:
#	rm -rf *.o cltest
#
#run:
#	adb push cltest /data/local/tmp; adb shell /data/local/tmp/cltest
