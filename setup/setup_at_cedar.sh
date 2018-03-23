#!/bin/bash

#
# T2KANA oscillation analysis code setup script for Compute Canada's cedar cluster (uses setup_at_cvmfs.sh)
#
# Input arguments: analysis code directory, framework directory
#

export T2KOSCANA="$(readlink -m $1)"
export VALOR_FMWK="$(readlink -m $2)"

source "${VALOR_FMWK}/setup/setup_at_cvmfs.sh" "${VALOR_FMWK}"
export LD_LIBRARY_PATH="${T2KOSCANA}/lib:${LD_LIBRARY_PATH}:${T2KOSCANA}/setup/lib_cedar"
export PATH="${T2KOSCANA}/bin:${T2KOSCANA}/production:${PATH}"

if [ ! -d "${T2KOSCANA}/setup/lib_cedar" ]; then
	mkdir -p $T2KOSCANA/setup/lib_cedar
	cd $T2KOSCANA/setup/lib_cedar
	ln -s /usr/lib64/libcom_err.so.2 libcom_err.so.2 
	ln -s /usr/lib64/libcrypto.so.10 libcrypto.so.10 
	ln -s /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib/libcrypto.so libcrypto.so.10.bak 
	ln -s /usr/lib64/libgssapi_krb5.so.2 libgssapi_krb5.so.2
	ln -s /usr/lib64/libk5crypto.so.3 libk5crypto.so.3 
	ln -s /usr/lib64/libkeyutils.so.1 libkeyutils.so.1 	
	ln -s /usr/lib64/libkrb5.so.3 libkrb5.so.3 
	ln -s /usr/lib64/libkrb5support.so.0 libkrb5support.so.0
	ln -s /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib/liblzma.so.5 liblzma.so.5 
	ln -s /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib/libpcre.so.1 libpcre.so.1
	ln -s /usr/lib64/libselinux.so.1 libselinux.so.1 
	ln -s /usr/lib64/libssl.so.10 libssl.so.10 
	ln -s /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib/libssl.so libssl.so.10.bak
	ln -s /cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/lib/libz.so.1 libz.so.1 
	cd -
fi

