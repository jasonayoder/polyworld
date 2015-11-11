conf=../../Makefile.conf
include ${conf}

target=${LIBRARY_TARGET}
blddir=${LIBRARY_BLDDIR}

cxxflags=-fPIC -I./ ${CXXFLAGS} ${OPENGL_CXXFLAGS} ${OMP_CXXFLAGS} ${CPPPROPS_CXXFLAGS}
ldflags=${SHARED_LDFLAGS}
libs=

include ${TARGET_MAK}