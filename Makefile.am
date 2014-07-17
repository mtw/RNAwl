bin_PROGRAMS = RNAwl
RNAwl_SOURCES =	main.c\
			moves.c\
			wl_options.c\
			wanglandau.c\
			wl_rna.c\
			wl_cmdline.c

AM_CFLAGS = ${GSL_CFLAGS} ${ViennaRNA_CFLAGS} -g3 -O0
AM_CPPFLAGS = -I${includedir} -I.

LDADD = ${GSL_LIBS} ${ViennaRNA_LIBS}
