####### Library info

LIBCPP = \
  latticelm.cc \
  data-lattice.cc \
  pylm.cc \
  sentence.cc \
  timer.cc \
  hierarchical-lm.cc \
  macros.cc \
  lexical-tm.cc \
  sampgen.cc


# AM_CXXFLAGS = $(BOOST_CPPFLAGS) $(OPENFST_CPPFLAGS) -I$(srcdir)/..
AM_CXXFLAGS = $(BOOST_CPPFLAGS) -I$(srcdir)/..  -I/home/oadams/tools/openfst-1.6.2/src/include

lib_LTLIBRARIES = liblatticelm.la

liblatticelm_la_SOURCES = ${LIBCPP}
liblatticelm_la_LDFLAGS = -version-info 0:0:0

####### Binary info

LDADD=liblatticelm.la \
    -L/home/oadams/tools/openfst-1.6.2/src/lib -lfst \
    $(BOOST_LDFLAGS) \
    $(BOOST_SYSTEM_LIB) \
    $(BOOST_PROGRAM_OPTIONS_LIB) \
    $(BOOST_SERIALIZATION_LIB) \
    $(BOOST_IOSTREAMS_LIB)

bin_PROGRAMS = latticelm lat1finalstate latsample

latticelm_SOURCES = latticelm-main.cc
latticelm_LDADD = $(LDADD)

lat1finalstate_SOURCES = lat1finalstate.cc
latsample_SOURCES = latsample.cc
