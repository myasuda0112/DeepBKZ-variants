CC	= g++-9
SRC	= main.cpp
CFLAGS	= -std=c++17 -march=native -mtune=native -O3
LDFLAGS	= -lgmp -lmpfr -lntl

all:
	${CC} ${CFLAGS} ${SRC} ${LDFLAGS}


