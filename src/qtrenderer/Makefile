conf=../../Makefile.conf
include ${conf}

target=${QTRENDERER_TARGET}
target_name=${QTRENDERER_TARGET_NAME}
blddir=${QTRENDERER_BLDDIR}

# The generated Qt project and make files
qt_project=./.qt.pro
qt_makefile=./.qt.mak

qt_modules=gui opengl core
cxxflags=${CXXFLAGS} ${LIBRARY_CXXFLAGS} ${OMP_CXXFLAGS} ${OPENGL_CXXFLAGS}
libs=${PWLIB_LDFLAGS} ${LIBRARY_LIBS} ${OMP_LIBS} ${OPENGL_LIBS}

.PHONY: ${qt_project} ${qt_makefile} qt_build clean

# Run the makefile generated by qmake.
qt_build: ${qt_makefile}
	+ make -f ${qt_makefile}
	@rm ${qt_makefile}
	@rm ${qt_project}

# Generate a makefile from the Qt project file
${qt_makefile}: ${qt_project}
	${QMAKE} -o ${qt_makefile} ${qt_project}
	$(call FIXUP_QT_MAKEFILE,$@)

# Generate the Qt project file
${qt_project}:
	@rm -f $@
	@echo "TEMPLATE = lib" >> $@
	@echo "CONFIG += qt" >> $@
	@echo "QT = ${qt_modules}" >> $@
	@echo "QMAKE_CXX=${CXX}" >> $@
	@echo "QMAKE_LINK=${LD}" >> $@
	@echo "QMAKE_CXXFLAGS += ${cxxflags}" >> $@
	@echo "QMAKE_LFLAGS_SONAME = ${QMAKE_LFLAGS_SONAME}" >> $@
	@echo "DEPENDPATH += ${LIBRARY_SRC}" >> $@
	@echo "LIBS += ${libs}" >> $@
	@echo "SOURCES += $(shell find . -name '*.cc')" >> $@
	@echo "HEADERS += $(shell find . -name '*.h')" >> $@
	@echo "TARGET = ${target_name}" >> $@
	@echo "DESTDIR = $(shell dirname ${target})" >> $@
	@echo "OBJECTS_DIR = ${blddir}" >> $@
	@echo "MOC_DIR = ${blddir}" >> $@

clean:
	rm -rf ${blddir}
	rm -f ${qt_makefile}
	rm -f ${qt_project}