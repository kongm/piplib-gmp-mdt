#############################################################################
SUBDIRS			=

#############################################################################
MAINTAINERCLEANFILES	= Makefile.in

INCLUDES		= -I$(top_builddir) -I$(top_srcdir)	\
			  -I$(top_builddir)/include	\
		          -I$(top_srcdir)/include

#############################################################################

lib_LTLIBRARIES		=
if WANT_LIBS
lib_LTLIBRARIES += libpiplib-gmp.la
libpiplib_gmp_la_SOURCES 	=			\
	integrer.c				\
	piplib.c				\
	tab.c					\
	sol.c					\
	traiter.c
libpiplib_gmp_la_CPPFLAGS = -DPIPMP_LINEAR_VALUE_IS_MP
endif

bin_PROGRAMS		=
if WANT_PIP
bin_PROGRAMS += pip-gmp
pip_gmp_SOURCES	= maind.c
pip_gmp_CPPFLAGS	= -DPIPMP_LINEAR_VALUE_IS_MP
pip_gmp_LDADD	= $(top_srcdir)/source/libpiplib-gmp.la
endif
AM_CFLAGS = -lm

