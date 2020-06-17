#
# Main VALOR/DUNE project Makefile.
#
# Author: Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(VALOR_EXPT)/make/Make.include

# define composite targets
#
BUILD_TARGETS =  print-make-info      \
		 make-bin-lib-dir     \
		 copy-scripts-to-bin  \
		 build-libs           \
		 build-apps

# define targets

all:     $(BUILD_TARGETS)
install:

print-make-info: FORCE
	@echo " "
	@echo " "
	@echo "***** Building VALOR/DUNE from source tree at: $(VALOR_EXPT_TOPDIR)"
	@echo "***** The source tree corresponds to VALOR version $(VALOR_VERSION)"
	@echo " "

copy-scripts-to-bin: | make-bin-lib-dir FORCE
	@echo " "
	@echo "** Copying scripts to bin..."
	cp ${VALOR_EXPT_TOPDIR}/make/valor_expt-config ${VALOR_EXPT_TOPDIR}/bin
	cp ${VALOR_EXPT_TOPDIR}/root/valor             ${VALOR_EXPT_TOPDIR}/bin
	cp ${VALOR_EXPT_TOPDIR}/root/pyvalor           ${VALOR_EXPT_TOPDIR}/bin

build-libs: make-bin-lib-dir copy-scripts-to-bin 

build-apps: make-bin-lib-dir build-libs copy-scripts-to-bin FORCE
	@echo " "
	@echo "** Building apps..."
	$(MAKE) -C ${VALOR_EXPT_TOPDIR}/src/dune/Apps;

make-bin-lib-dir: FORCE
	@echo " "
	@echo "** Creating VALOR/DUNE lib and bin directories..."
	cd ${VALOR_EXPT_TOPDIR}; \
	[ -d bin ] || mkdir bin; chmod 755 bin;\
	[ -d lib ] || mkdir lib; chmod 755 lib;

# copy-rootdict-to-lib: make-bin-lib-dir build-libs FORCE
# 	@echo " "
# 	@echo "** Copying ROOT dictionary PCM files to ${VALOR_EXPT_TOPDIR}/lib"
# 	$(shell find "${VALOR_EXPT_TOPDIR}/src/dune/" -name "_RootDict*_rdict.pcm" | xargs -n 1 -I{} cp {} "${VALOR_EXPT_TOPDIR}/lib/" )

purge: FORCE
	@echo " "
	@echo "** Purging..."
	$(MAKE) -C src/dune/Apps    purge 

clean: clean-files clean-dir clean-etc

force-clean: clean-files force-clean-dir clean-etc

clean-files: FORCE
	@echo " "
	@echo "** Cleaning..."
	rm -f ${VALOR_EXPT_TOPDIR}/bin/valor_expt-config; \
	rm -f ${VALOR_EXPT_TOPDIR}/bin/pyvalor; \
	rm -f ${VALOR_EXPT_TOPDIR}/bin/valor; \
	$(MAKE) -C src/dune/Apps    clean; \
	rm -f ${VALOR_EXPT_TOPDIR}/src/dune/Conventions/Build.h; \
	rm -f ${VALOR_EXPT_TOPDIR}/src/dune/Conventions/Version.h; \
	cd $(VALOR_EXPT_TOPDIR);

clean-dir: FORCE
	@echo "Deleting VALOR/DUNE lib and bin directories...";\
	cd $(VALOR_EXPT_TOPDIR);\
	[ ! -d ./bin ] || rmdir ./bin;\
	[ ! -d ./lib ] || rmdir ./lib;\
	[ ! -d ./src/dune/Conventions ] || rmdir ./src/dune/Conventions

force-clean-dir: FORCE
	@echo "Deleting VALOR/DUNE lib and bin directories...";
	cd $(VALOR_EXPT_TOPDIR);
	rm -rf ./bin;
	rm -rf ./lib;
	rm -rf ./src/dune/Conventions

clean-etc: FORCE
	cd $(VALOR_EXPT_TOPDIR); \
	rm -f ./*log; \
	cd ${VALOR_EXPT_TOPDIR}

FORCE:
