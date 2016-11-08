SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

include $(T2KOSCANA)/make/Make.include

# define composite targets 
#
BUILD_TARGETS =    print-make-info     \
                   make-bin-lib-dir    \
                   copy-scripts-to-bin \
                   build-fitter        \
                   build-fit-apps      \
                   build-test-apps     

# define targets
all:     $(BUILD_TARGETS)


print-make-info: FORCE
	@echo " "
	@echo " "
	@echo "***** Building VALOR/T2XK from source tree at: $(T2KOSCANA)"
	@echo "***** The source tree corresponds to VALOR version $(VALOR_VERSION)"
	@echo " "

copy-scripts-to-bin: make-bin-lib-dir FORCE
	@echo " "
	@echo "** Copying scripts to bin..."
	cp ${T2KOSCANA}/make/valor_t2xk-config ${T2KOSCANA}/bin
	cp ${T2KOSCANA}/root/root2k ${T2KOSCANA}/bin
	cp ${T2KOSCANA}/root/pyroot2k ${T2KOSCANA}/bin

build-fitter: make-bin-lib-dir FORCE
	@echo " "
	@echo "** Building T2K oscillation fitter library..."
	$(MAKE) -C ${T2KOSCANA}/fitter

make-bin-lib-dir: FORCE
	@echo " "
	@echo "** Creating VALOR/T2XK lib and bin directories..."
	cd ${T2KOSCANA}; \
	[ -d bin ] || mkdir bin; chmod 755 bin;\
	[ -d lib ] || mkdir lib; chmod 755 lib;

build-fit-apps: make-bin-lib-dir build-fitter FORCE
	@echo " "
	@echo "** Building T2XK fit apps..."
	$(MAKE) -C ${T2KOSCANA}/apps/fit/

build-test-apps: make-bin-lib-dir build-fitter FORCE
	@echo " "
	@echo "** Building T2XK test apps..."
	$(MAKE) -C ${T2KOSCANA}/apps/test/

clean: clean-files clean-dir

force-clean: clean-files force-clean-dir

clean-files: FORCE
	@echo " "
	@echo "** Cleaning-up..."
	rm -f ${T2KOSCANA}/bin/valor_t2xk-config; \
	rm -f ${T2KOSCANA}/bin/pyroot2k; \
	rm -f ${T2KOSCANA}/bin/root2k; \
	make clean -C ${T2KOSCANA}/fitter;\
	make clean -C ${T2KOSCANA}/apps/fit/;\
	make clean -C ${T2KOSCANA}/apps/test/;
	@echo " "
	@echo " "
	@echo " "
	@echo "Deleting libaries in $T2KOSCANA/fitter. They have moved, so this will be disabled in the near future"
	rm -f ${T2KOSCANA}/fitter/libT2KNuOsc*
	@echo " "
	@echo " "
	@echo " "


clean-dir: FORCE
	@echo "Deleting VALOR/T2XK lib and bin directories...";\
	cd $(T2KOSCANA);\
	[ ! -d ./bin ] || rmdir ./bin;\
	[ ! -d ./lib ] || rmdir ./lib 

force-clean-dir: FORCE
	@echo "Deleting VALOR/T2XK lib and bin directories...";\
	cd $(VALOR_LBNX);\
	[ ! -d ./bin ] || rm -rf ./bin;\
	[ ! -d ./lib ] || rm -rf ./lib 

FORCE:

