#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging

"""
This script will run the macro test_systematic_error_propagation.C for all systematic possiblities, producing
a latex file (and then pdf) with the outputs. Things you can change are:
 - conffile - configuration file
 - basedir - directory to save tex & pdf. Default is 'current working directory'
 - filename - output filename. Default is systematic_tables.tex
 - systematics - the systematics you want to run over (defaults to all currently known)
 - samples - the SK samples you want to consider (default 1Rmu)
 - oscparam - the oscillation parameters you want to calculate shifts relative to (defaults to global best fits - you can specify more than 1 set of parameters to run over)
            - (Note that if sinsqtheta13 is identical in 2 cases, you'll need to add another distinguisher to the label definition at the end of root macro, otherwise the tex won't compile)
 - spectralreco - include ccqe ereco spectral information? (default true - you can specify both true and false)
 - spectraltrue - include ccqe etrue spectral information? (default true - you can specify both true and false)
 - banffbest - run with banff bestfit? (default true - you can specify both true and false)

Run this by /path/to/test_systematic_error_propagation_runAll.py (you may need to chmod +x)
It will then create the .tex, .pdf, and associated files in the current working directory
"""

conffile = "analysis_configuration.xml"

#will create directories relative to the current working directory, if they don't already exist
basedir = os.getenv("PWD")
filename = basedir+"/systematic_tables.tex"
if os.path.isfile(filename):
   sys.exit("systematic_tables.tex already exists. Will not overwrite")
f = open(filename,'w')
f.write("\\documentclass[10pt,a4paper]{article}\n")
f.write("\\usepackage{amsmath}\n")
f.write("\\usepackage{graphicx}\n")
f.write("\\usepackage{dcolumn}\n")
f.write("\\begin{document}\n")
f.close()

#specify the cuts you want to apply here
systematics = []
#systematics.append("sk_1Rmu_ccqe_ring_counting_eff") #not supported in the macro
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_0")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_1")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_2")
systematics.append("sk_1Rmu_ccqe_other_eff")
systematics.append("sk_1Rmu_ccnonqe_eff")
systematics.append("sk_1Rmu_nc_eff")
systematics.append("sk_1Rmu_nuecc_eff")
systematics.append("sk_energy_scale")
systematics.append("banff") #not supported in this macro
"""
for i in xrange(26):
   systematics.append("banff_" + '%(number)02d' % {"number": i})
"""

samples = []
samples.append("1Rmu")
#samples.append("1Re")
   
oscparams = [] #fill with oscillation parameters, in order mixing12, theta13, theta23, deltacp, dmsq12, dmsq23hypothesis 
oscparams.append([0.86, 0.1,  1,  0,  7.58E-5,  2.32E-3])
#oscparams.append([0.86, 0,    1,  0,  7.58E-5,  2.32E-3])

spectralreco = []
spectralreco.append("true")
#spectralreco.append("false")

spectraltrue = []
spectraltrue.append("true")
#spectraltrue.append("false")

banffbest = []
banffbest.append("true")
banffbest.append("false")

for osc in oscparams:

   for spr in spectralreco:

      for spt in spectraltrue:

         for sam in samples:
            
            for sys in systematics:
          
               for ban in banffbest:
                  
                  command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_systematic_error_propagation.C("'+sys+'","'+sam+'","3f_osc",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+','+ban+',true,'+spr+','+spt+',"'+filename+'","'+conffile+'")\''
                  print command
                  os.system(command)
          

f = open(filename,'a')
f.write("\\end{document}\n")
f.close()

command = 'pdflatex '+filename
print command
os.system(command)
