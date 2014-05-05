bin_PROGRAMS = wanglandau
wanglandau_SOURCES =	main.c\
			moves.c\
			wl_options.c\
			wanglandau.c\
			wl_rna.c

AM_CFLAGS = -g3 -O0 ${GSL_CFLAGS} ${ViennaRNA_CFLAGS}
INCLUDES= -I${includedir} -I.

LDADD = ${GSL_LIBS} ${ViennaRNA_LIBS}
