EXEC_NAME=Main
CC=g++

CFLAGS = -O3 -std=c++17
#FLAGS = -Wall -Wextra -Wno-unused-result -Wno-char-subscripts -Wshadow -Wfloat-equal -Wconversion -Wformat-signedness -Wvla -Wduplicated-cond -Wlogical-op -Wredundant-decls -ggdb3 -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2 -fsanitize=undefined,address,float-divide-by-zero,float-cast-overflow -fno-omit-frame-pointer -fno-optimize-sibling-calls -fstack-protector-all -fno-sanitize-recover=all -O2 

SRC =                         \
	nauty/nauty.h	              \
	nauty/nauty.c	              \
	nauty/nautil.c	            \
	nauty/naugraph.c	     	    \
	nauty/schreier.c	          \
	nauty/naurng.c              \
	ESU.cpp                     \
	GTrie.cpp                   \
	Hypergraph.cpp              \
	Isomorphism.cpp             \
	Main.cpp

OBJ =  ${SRC:.cpp=.o}

#------------------------------------------------------------

all: ${EXEC_NAME}

${EXEC_NAME}: ${OBJ}
	${CC} ${CFLAGS} ${CLIBS} -o ${EXEC_NAME} ${OBJ}

%.o: %.cpp
	${CC} ${CFLAGS} -c -o $@ $+

clean:
	rm ${EXEC_NAME} *.o *~ *# -rf




