SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

include $(VALOR_FMWK)/make/Make.include

# define composite targets
#
BUILD_TARGETS =    print-make-info     \
                   make-bin-lib-dir    \
                   copy-scripts-to-bin \
                   build-fitter        \
                   build-fit-apps

# define targets
all:     $(BUILD_TARGETS)


print-make-info: FORCE
	@echo " "
	@echo " "
	@echo "***** Building VALOR/T2XK from source tree at: $(T2KOSCANA)"
	@echo "***** The source tree corresponds to VALOR version $(VALOR_VERSION)"
	@echo " "

copy-scripts-to-bin: FORCE
	@echo " "
	@echo "** Copying scripts to bin..."
	cp ${T2KOSCANA}/make/valor_t2xk-config ${T2KOSCANA}/bin

build-fitter: FORCE
	@echo " "
	@echo "** Building T2K oscillation fitter library..."
	cd ${T2KOSCANA}/fitter;\
	make;\
	cd ${T2KOSCANA};

make-bin-lib-dir: FORCE
	@echo " "
	@echo "** Creating VALOR/T2XK lib and bin directories..."
	cd ${T2KOSCANA}; \
	[ -d bin ] || mkdir bin; chmod 755 bin;\
	[ -d lib ] || mkdir lib; chmod 755 lib;

build-fit-apps: FORCE
	@echo " "
	@echo "** Building T2XK fit apps..."
	cd ${T2KOSCANA}/apps/fit/;\
	make;\
	cd ${T2KOSCANA};


clean: clean-files clean-dir

force-clean: clean-files force-clean-dir

clean-files: FORCE
	@echo " "
	@echo "** Cleaning-up..."
	rm -f ${T2KOSCANA}/bin/valor_t2xk-config; \
	cd ${T2KOSCANA}/fitter;\
	make clean;\
	cd ${T2KOSCANA}/apps/fit/;\
	make clean;\
	cd ${T2KOSCANA};

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

