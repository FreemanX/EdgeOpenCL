all: tinymembench

CC=aarch64-linux-android-g++
#CC=aarch64-linux-android-clang++
CFLAGS= -I./include -L./lib -g -O1 -fPIC -fPIE -static-libstdc++ -pie -fuse-ld=gold 

tinymembench: main.c util.o util.h asm-opt.h version.h asm-opt.o x86-sse2.o arm-neon.o mips-32.o aarch64-asm.o
	${CC} ${CFLAGS} -o tinymembench main.c util.o asm-opt.o x86-sse2.o arm-neon.o mips-32.o aarch64-asm.o -lm

util.o: util.c util.h
	${CC} ${CFLAGS} -c util.c

asm-opt.o: asm-opt.c asm-opt.h x86-sse2.h arm-neon.h mips-32.h
	${CC} ${CFLAGS} -c asm-opt.c

x86-sse2.o: x86-sse2.S
	${CC} ${CFLAGS} -c x86-sse2.S

arm-neon.o: arm-neon.S
	${CC} ${CFLAGS} -c arm-neon.S

aarch64-asm.o: aarch64-asm.S
	${CC} ${CFLAGS} -c aarch64-asm.S

mips-32.o: mips-32.S
	${CC} ${CFLAGS} -c mips-32.S

clean:
	-rm -f tinymembench
	-rm -f tinymembench.exe
	-rm -f *.o
