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
 - filename - output filename (defaults to 'current working directory'/systematic_tables.tex
 - systematics - the systematics you want to run over (defaults to all currently known)
 - samples - the SK samples you want to consider (default 1RmuT, 1ReT)
 - oscparam - the oscillation parameters you want to calculate shifts relative to (defaults to global best fits - you can specify more than 1 set of parameters to run over)
            - (Note that if sinsqtheta13 is identical in 2 cases, you'll need to add another distinguisher to the label definition at the end of root macro, otherwise the tex won't compile)
 - spectrals - include ccqe spectral information (default true - you can specify both true and false)

Run this by /path/to/test_systematic_error_propagation_runAll.py (you may need to chmod +x)
It will then create the .tex, .pdf, and associated files in the current working directory
"""

#will create directories relative to the current working directory, if they don't already exist
base = os.getenv("PWD")
filename = base+"/systematic_tables.tex"
if os.path.isfile(filename):
   sys.exit("systematic_tables.tex already exists. Will not overwrite")
f = open(filename,'w')
f.write("\\documentclass[10pt,a4paper]{article}\n")
f.write("\\usepackage{amsmath}\n")
f.write("\\usepackage{graphicx}\n")
f.write("\\usepackage{dcolumn}\n")
#f.write("\n")
f.write("\\begin{document}\n")
f.close()

#specify the cuts you want to apply here
systematics = []
systematics.append("sk_energy_scale")
systematics.append("sk_ccqe_other_eff")
#systematics.append("sk_ccqe_ring_counting_eff") #not supported in the macro
systematics.append("sk_ccqe_ring_counting_eff_0")
systematics.append("sk_ccqe_ring_counting_eff_1")
systematics.append("sk_ccqe_ring_counting_eff_2")
systematics.append("sk_ccnonqe_eff")
systematics.append("sk_nc_eff")
systematics.append("sk_nuecc_eff")
systematics.append("sk_1Resig_eff")
systematics.append("sk_1Rebkg_eff")
systematics.append("nd280_eff")
systematics.append("nd280_norm")
systematics.append("flux_norm")
#systematics.append("flux_shape") #not supported in the macro
systematics.append("xsec_ccqe")
systematics.append("xsec_cc1pi_to_ccqe")
systematics.append("xsec_cccoh_to_ccqe")
systematics.append("xsec_ccoth_to_ccqe")
systematics.append("xsec_nc_to_ccqe")
systematics.append("xsec_nc1pi0_to_ccqe")
systematics.append("xsec_nccoh_to_ccqe")
systematics.append("xsec_ncoth_to_ccqe")
systematics.append("xsec_nue_to_numu")
systematics.append("nuclear_effect")

samples = []
samples.append("1Rmu")
samples.append("1Re")
   
oscparams = [] #fill with oscillation parameters, in order mixing12, theta13, theta23, deltacp, dmsq12, dmsq23hypothesis 
oscparams.append([0.86, 0.1,  1,  0,  7.58E-5,  2.32E-3])
#oscparams.append([0.86, 0,    1,  0,  7.58E-5,  2.32E-3])

spectrals = []
spectrals.append("true")
#spectrals.append("false")


for osc in oscparams:

   for spe in spectrals:

    for sam in samples:

       for sys in systematics:

          command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_systematic_error_propagation_nue.C("'+sys+'","'+sam+'","3f_osc",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+',true,'+spe+',"'+filename+'")\''
          print command
          os.system(command)
          

f = open(filename,'a')
f.write("\\end{document}\n")
f.close()

command = 'pdflatex '+filename
print command
os.system(command)
