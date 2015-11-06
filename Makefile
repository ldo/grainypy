# This Makefile is now used for testing only.

PYTHONVER=python3.4m
CFLAGS=-g -I/usr/include/${PYTHONVER} -fPIC -Wall

grainyx.so : grainyx.o
	$(CC) $^ -L${PYTHONVER}/config -l${PYTHONVER} -lpng -shared -o $@

grainyx.o : grainyx.c

clean :
	rm -f grainyx.so grainyx.o build/

.PHONY : clean
