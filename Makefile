SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

all: build-fitter

build-fitter: FORCE
	@echo " "
	@echo "** Building T2K oscillation fitter library..."
	cd fitter; make; cd ..

clean: FORCE
	@echo " "
	@echo "** Cleaning-up..."
	cd fitter; make clean; cd ..

FORCE:

