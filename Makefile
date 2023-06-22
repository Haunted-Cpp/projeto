EXEC_NAME=hypermotif
CC=g++

CFLAGS = -O3 -std=c++17 -w
#CFLAGS = -std=c++17 -Wall -Wextra -Wno-unused-result -Wno-char-subscripts -Wshadow -Wfloat-equal -Wconversion -Wformat-signedness -Wvla -Wduplicated-cond -Wlogical-op -Wredundant-decls -ggdb3 -fno-optimize-sibling-calls -fstack-protector-all -fno-sanitize-recover=all -O3 -fsanitize=undefined,address,float-divide-by-zero,float-cast-overflow -fno-omit-frame-pointer -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2  -D_GLIBCXX_DEBUG 

SRC =                         \
	nauty/nauty.h	              \
	nauty/nauty.c	              \
	nauty/nautil.c	            \
	nauty/naugraph.c	     	    \
	nauty/schreier.c	          \
	nauty/naurng.c              \
	ESU.cpp                     \
	Hypergraph.cpp              \
	IsomorphismHyper.cpp             \
	FaSE/Fase.cpp 	\
	FaSE/Label.cpp	\
	FaSE/IGtrie.cpp	\
	FaSE/Isomorphism.cpp	\
	FaSE/Timer.cpp	\
	FaSE/DynamicGraph.cpp\
	FaSE/GraphMatrix.cpp	\
	FaSE/GraphUtils.cpp	\
	FaSE/Random.cpp	\
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
	cd FaSE && rm ${EXEC_NAME} *.o *~ *# -rf



