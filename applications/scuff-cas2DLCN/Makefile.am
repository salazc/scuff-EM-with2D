#SUBDIRS = libTDRT libAmosBessel

bin_PROGRAMS = scuff-cas2DLCN

scuff_cas2DLCN_SOURCES = 		\
 scuff-cas2DLCN.h  		\
 scuff-cas2DLCN.cc			\
 CreateC2DWorkspace.cc		\
 CubatureMethods.cc		\
 ImageObjects.cc		\
 XQIntegrand.cc

scuff_cas2DLCN_LDADD = \
 $(top_builddir)/libs/libTDRT/libTDRT.la				\
 $(top_builddir)/libs/libAmosBessel/libAmosBessel.la		\
 $(top_builddir)/libs/libSGJC/libSGJC.la		\
 $(top_builddir)/libs/libMDInterp/libMDInterp.la	\
 $(top_builddir)/libs/libMatProp/libMatProp.la	\
 $(top_builddir)/libs/libhmat/libhmat.la		\
 $(top_builddir)/libs/libhrutil/libhrutil.la


AM_CPPFLAGS = \
  -I$(top_srcdir)/libs/libhrutil				\
  -I$(top_srcdir)/libs/libTDRT		\
  -I$(top_srcdir)/libs/libAmosBessel 	\
  -I$(top_srcdir)/libs/libSGJC       \
  -I$(top_srcdir)/libs/libMDInterp   \
  -I$(top_srcdir)/libs/libMatProp    \
  -I$(top_srcdir)/libs/libhmat
