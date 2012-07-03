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

mixing23s = [1.0]
dm23sqs = [0.00235]
dm12sq = -99
mixing13 = -99
mixing12 = -99
deltacp = -99
banffs = ["false","true"]


for mixing23 in mixing23s:
    for dm23sq in dm23sqs:
        
        for banff in banffs:
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_sk_ereco_spectrum_nominal.C("1Rmu","3f_osc",'+str(mixing12)+','+str(mixing13)+','+str(mixing23)+','+str(deltacp)+','+str(dm12sq)+','+str(dm23sq)+','+banff+',"analysis_configuration.xml",true)\''
            print command
            os.system(command)

        command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_systematic_error_propagation_all.C("1Rmu",'+str(mixing12)+','+str(mixing13)+','+str(mixing23)+','+str(deltacp)+','+str(dm12sq)+','+str(dm23sq)+',"analysis_configuration.xml")\''
        print command
        os.system(command)


# not oscillated
for banff in banffs:
    command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_sk_ereco_spectrum_nominal.C("1Rmu","3f_osc",0,0,0,0,0,0,'+banff+',"analysis_configuration.xml",true)\''
    print command
    os.system(command)

command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_systematic_error_propagation_all.C("1Rmu",0,0,0,0,0,0,"analysis_configuration.xml")\''
print command
os.system(command)

#expected 1 sigma deviations
command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_systematic_error_propagation_expectation.C(true,"systematic_tables_expectation.tex","analysis_configuration.xml")\''
print command
os.system(command)
