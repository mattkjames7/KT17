
all: o32 o64 libkt17_32.so libkt17_64.so rem_obj
	
libkt17_32.so:
	gcc -m32 -shared -o libkt17_32.so kt17.o32 kt17trace.o32 -lm

libkt17_64.so:
	gcc -m64 -shared -o libkt17_64.so kt17.o64 kt17trace.o64 -lm
	 
o32:
	gcc -c -Wall -Werror -fpic -m32 kt17.cc -o kt17.o32
	gcc -c -Wall -Werror -fpic -m32 kt17trace.cc -o kt17trace.o32

o64:
	gcc -c -Wall -Werror -fpic -m64 kt17.cc -o kt17.o64
	gcc -c -Wall -Werror -fpic -m64 kt17trace.cc -o kt17trace.o64

rem_obj:
	rm -vf *.o*

clean:
	rm -vf *.so
