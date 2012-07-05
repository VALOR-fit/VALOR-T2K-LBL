#! /usr/bin/env python

"""
Run this to produce a number of 

"""

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging

optionParser = optparse.OptionParser(usage ="""\
usage: %prog [options]
    Run many root2k scripts for comparisons in the numu oscillation analysis""")
#names of things
optionParser.add_option("-n","--no_nsk", action="store_true", dest="skip_nsk", default=False, help="If set, will not run test_nsk.C")
optionParser.add_option("-s","--no_spect_nom", action="store_true", dest="skip_spect_nom", default=False, help="If set, will not run test_sk_ereco_spectrum_nominal.C")
optionParser.add_option("-e","--no_error_prop", action="store_true", dest="skip_error_prop", default=False, help="If set, will not run test_systematic_error_propagation_all.C")
(opt, args) = optionParser.parse_args()

mixing23s = [1.0]
dm23sqs = [0.00240]
dm12sq = -99
mixing13 = -99
mixing12 = -99
deltacp = -99
banffs = ["false","true"]

for mixing23 in mixing23s:
    for dm23sq in dm23sqs:
        
        for banff in banffs:
            #nominal spectra
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_sk_ereco_spectrum_nominal.C("1Rmu","3f_osc",'+str(mixing12)+','+str(mixing13)+','+str(mixing23)+','+str(deltacp)+','+str(dm12sq)+','+str(dm23sq)+','+banff+',"analysis_configuration.xml",true)\''
            print command
            if not opt.skip_spect_nom:
                os.system(command)

        #1 sigma errors
        command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_systematic_error_propagation_all.C("1Rmu",'+str(mixing12)+','+str(mixing13)+','+str(mixing23)+','+str(deltacp)+','+str(dm12sq)+','+str(dm23sq)+',"analysis_configuration.xml")\''
        print command
        if not opt.skip_error_prop:
            os.system(command)

        #nsk (table summarising banff/no banff & osc/no osc)
        command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_nsk.C\''+str(mixing12)+','+str(mixing13)+','+str(mixing23)+','+str(deltacp)+','+str(dm12sq)+','+str(dm23sq)+'"1Rmu","analysis_configuration.xml")\''
        print command
        if not opt.skip_nsk:
            os.system(command)

# not oscillated #############################################

#nominal spectra
for banff in banffs:
    command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_sk_ereco_spectrum_nominal.C("1Rmu","3f_osc",0,0,0,0,0,0,'+banff+',"analysis_configuration.xml",true)\''
    print command
    if not opt.skip_spect_nom:
        os.system(command)

#1 sigma errors
command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_systematic_error_propagation_all.C("1Rmu",0,0,0,0,0,0,"analysis_configuration.xml")\''
print command
if not opt.skip_error_prop:
    os.system(command)




#expected 1 sigma deviations
command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_systematic_error_propagation_expectation.C(true,"systematic_tables_expectation.tex","analysis_configuration.xml")\''
print command
if not opt.skip_error_prop:
    os.system(command)
