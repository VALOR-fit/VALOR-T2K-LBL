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
(opt, args) = optionParser.parse_args()

conffile = "analysis_configuration.xml"

save_toy_expts = "false"

modes1 = [0,1,2,3,4] #combinations
modes3 = range(200,215) #unconstranied xsec
modes4 = range(300,306) #sk detector

modes = modes1
if opt.do_all:
    modes = modes + modes3 + modes4

oscparams = [] #fill with oscillation parameters, in order mixing12, theta13, theta23, deltacp, dmsq12, dmsq23hypothesis
oscparams.append([  0,   0,   0,   0,   0,   0]) #no oscillations!
oscparams.append([-99, -99, -99, -99, -99, -99]) #use global best fits!
#oscparams.append([0.86, 0.1,  1,  0,  7.58E-5,  2.32E-3])
#oscparams.append([0.86, 0,    1,  0,  7.58E-5,  2.32E-3])


for osc in oscparams:

    update = "false"
    
    for mode in modes:

        command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_sk_ereco_spectrum_twk_sys_throw.C('+str(mode)+','+update+',"3f_osc",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+','+save_toy_expts+',"'+conffile+'")\''
        print command
        if os.system(command) == 1:
            sys.exit(1)
        if update == "false":
            update = "true"
