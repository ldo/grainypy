# This Makefile is now used for testing only.

PYTHONVER=python3.4m
CFLAGS=-g -I/usr/include/${PYTHONVER} -fPIC -Wall -Wno-parentheses

grainyx.so : grainyx.o
	$(CC) $^ -L${PYTHONVER}/config -l${PYTHONVER} -lcairo -shared -o $@

grainyx.o : grainyx.c

clean :
	rm -rf grainyx.so grainyx.o build/

.PHONY : clean
