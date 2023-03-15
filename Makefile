COMPILER = g++
# FLAGS = -Wall -Wextra -Wno-unused-result -Wno-char-subscripts -Wshadow -Wfloat-equal -Wconversion -Wformat-signedness -Wvla -Wduplicated-cond -Wlogical-op -Wredundant-decls -ggdb3 -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2 -fsanitize=undefined,address,float-divide-by-zero,float-cast-overflow -fno-omit-frame-pointer -fno-optimize-sibling-calls -fstack-protector-all -fno-sanitize-recover=all -O2 
#FLAGS = -Wall -Wextra -Wshadow -Ofast

FLAGS = -Ofast

VERSION = -std=c++17

NAUTY_SRC =                   \
	nauty/nauty.h	              \
	nauty/nauty.c	              \
	nauty/nautil.c	            \
	nauty/naugraph.c	     	    \
	nauty/schreier.c	          \
	nauty/naurng.c

output: Main.o Isomorphism.o Hypergraph.o
	$(COMPILER) -o main Main.o Isomorphism.o Hypergraph.o $(NAUTY_SRC) $(FLAGS) $(VERSION)

Main.o: Main.cpp
	$(COMPILER) -c Main.cpp $(NAUTY_SRC) $(FLAGS) $(VERSION)

Isomorphism.o: Isomorphism.cpp Isomorphism.hpp
	$(COMPILER) -c Isomorphism.cpp $(NAUTY_SRC) $(FLAGS) $(VERSION)

Hypergraph.o: Hypergraph.cpp Hypergraph.hpp
	$(COMPILER) -c Hypergraph.cpp $(NAUTY_SRC) $(FLAGS) $(VERSION)

clean:
	rm -f *.o main
