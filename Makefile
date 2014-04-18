CC = gcc
VRNA_INC = $(shell pkg-config --cflags "RNAlib2 >= 2.2" gsl)
VRNA_LIB = $(shell pkg-config --libs "RNAlib2 >= 2.2" gsl)

#GSL_INC = $(shell pkg-config --cflags gsl)
#GSL_LIB = $(shell pkg-config --libs gsl)

INCLUDES =  $(VRNA_INC) -DLOOP_EN
LIBS    =  ${VRNA_LIB} -lm
CFLAGS = -Wall -g3
# rules

OFILES =        moves.o
EXEFILE =       moves

all:                            $(EXEFILE)

debug:				CFLAGS := $(CFLAGS) -g3

debug:				all

debug_parallel:			CFLAGS := $(CFLAGS) -g3 -DUSE_OPENMP -fopenmp

debug_parallel:			all

profile:			CFLAGS := $(CFLAGS) -pg

profile:			all

profile_parallel:		CFLAGS := $(CFLAGS) -g3 -pg -DUSE_OPENMP -fopenmp

profile_parallel:		all

parallel:			CFLAGS := $(CFLAGS) -DUSE_OPENMP -fopenmp

parallel:			all

all:                            $(EXEFILE)

clean:				
				rm -f $(OFILES)

moves.o:			moves.c moves.h
				$(CC)  -c $(CFLAGS) $(INCLUDES) $(GSL_INC) $(FOO)  moves.c

moves:	 			$(OFILES)
				$(CC) $(CFLAGS) -o $(EXEFILE) $(OFILES) $(LIBS)

# End of file
