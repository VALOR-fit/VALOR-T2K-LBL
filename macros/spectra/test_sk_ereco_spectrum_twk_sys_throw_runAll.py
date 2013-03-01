#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging

"""
This script will run the macro test_sk_ereco_spectrum_twk_sys_throw.C for all systematic possiblities, producing
a ROOT file with the outputs. Things you can change are:
 - conffile - configuration file
 - oscparam - the oscillation parameters you want to calculate shifts relative to (defaults to global best fits - you can specify more than 1 set of parameters to run over)
            - (Note that if sinsqtheta13 is identical in 2 cases, you'll need to add another distinguisher to the label definition at the end of root macro, otherwise the tex won't compile)
 - save_toy_expts - save the toy experiment ereco spectra? Warning: will produce big files! (20000 toy spectra per syst!)
 
Run this by /path/to/test_sk_ereco_spectrum_twk_sys_throw_runAll.py [--all] (you may need to chmod +x)
[if the --all option is set, will run systematics individually also, on top of the normal 5 groups]
It will then create the .root in the current working directory
"""

optionParser = optparse.OptionParser(usage ="""\
usage: %prog [options]
    Run test_sk_ereco_spectrum_twk_sys_throw.C many times""")
#names of things
optionParser.add_option("-a","--all", action="store_true", dest="do_all", default=False, help="If set, will run indiviudal systematics, as well as the categories")
optionParser.add_option("-n","--nexp", action="store", dest="nexp", type="int", default=-1, help="How many experiments to generate at each point. Required")
optionParser.add_option("-s","--spectra", action="store_true", dest="save_spectra", default=False, help="If set, will save the spectra. Requires PRODUCTIONMODE to be off")
(opt, args) = optionParser.parse_args()

nexp = opt.nexp
if nexp <= 0:
    sys.exit("specify number of events")
    
conffile = "analysis_configuration.xml"

save_toy_expts = "false"

no_spectra = "true"
if opt.save_spectra:
    no_spectra = "false"

modes1 = [0,2,3,4,5] #combinations
modes2 = range(100,125) #banff
modes3 = range(200,216) #unconstranied xsec
modes4 = range(300,309) #sk detector
modes5 = range(400,406) #fsi
modes6 = [7,8] #banff flux / banff xsec individually

modes = modes1
if opt.do_all:
    modes = modes + modes2 + modes3 + modes4 + modes5 + modes6

oscparams = [] #fill with oscillation parameters, in order mixing12, theta13, theta23, deltacp, dmsq12, dmsq23hypothesis
oscparams.append([  0,   0,   0,   0,   0,   0]) #no oscillations!
oscparams.append([-99, -99, -99, -99, -99, -99]) #use global best fits!
#oscparams.append([0.86, 0.1,  1,  0,  7.58E-5,  2.32E-3])
#oscparams.append([0.86, 0,    1,  0,  7.58E-5,  2.32E-3])

for osc in oscparams:

    update = "false"
    
    for mode in modes:

        command = 'root2k -l -b -q \'$T2KOSCANA/macros/spectra/test_sk_ereco_spectrum_twk_sys_throw.C('+str(mode)+','+update+',"3f_osc",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+','+save_toy_expts+',"'+conffile+'",'+str(nexp)+','+no_spectra+')\''
        print command
        if os.system(command) == 1:
            sys.exit(1)
        if update == "false":
            update = "true"

        #don't forget the banff prefit
        if mode == 5 or mode == 0 or (mode >= 100 and mode < 200):
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/spectra/test_sk_ereco_spectrum_twk_sys_throw.C('+str(mode)+','+update+',"3f_osc",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+','+save_toy_expts+',"analysis_configuration_2012a_run1+2+3_neut_finalBANFF_coarseErecoSplines_6binFSI+SI_MAQE_prefit.xml",'+str(nexp)+','+no_spectra+')\''
            print command
            if os.system(command) == 1:
                sys.exit(1)

            
