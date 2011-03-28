SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

all: fitter

fitter: FORCE
	@echo " "
	@echo "** Building T2K oscillation fitter..."
	cd fitter; make; cd ..

clean: FORCE
	@echo " "
	@echo "** Cleaning-up..."
	cd fitter; make clean; cd ..

FORCE:

