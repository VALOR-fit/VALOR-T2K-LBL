#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging

"""
This script will run the macro test_sk_ereco_spectrum_twk_sys.C for all systematic possiblities, producing
a ROOT file with the outputs (plus lots of .png's). Things you can change are:
 - conffile - configuration file
 - oscparam - the oscillation parameters you want to calculate shifts relative to (defaults to global best fits - you can specify more than 1 set of parameters to run over)
            - (Note that if sinsqtheta13 is identical in 2 cases, you'll need to add another distinguisher to the label definition at the end of root macro, otherwise the tex won't compile)
 - emax     - maximum of the Ereco (x-axis) to plot in the images
 
Run this by /path/to/test_sk_ereco_spectrum_twk_sys_runAll.py (you may need to chmod +x)
It will then create the .root, *.png, in the current working directory
"""

conffile = "analysis_configuration.xml"
emax = 7

modes1 = [0,1,2,3,4,5,6,7,8,9,10] #combinations
modes2 = range(100,139) #banff
modes3 = range(200,215) #unconstranied xsec
modes4 = range(300,308) #sk detector
modes5 = range(400,879) #fsi
modes6 = range(500,504) #special comparison tweaks

modes = modes1 + modes2 + modes3 + modes4 + modes5 #+ modes6

oscparams = [] #fill with oscillation parameters, in order mixing12, theta13, theta23, deltacp, dmsq12, dmsq23hypothesis
oscparams.append([  0,   0,   0,   0,   0,   0]) #no oscillations!
oscparams.append([-99, -99, -99, -99, -99, -99]) #use global best fits!
#oscparams.append([0.86, 0.1,  1,  0,  7.58E-5,  2.32E-3])
#oscparams.append([0.86, 0,    1,  0,  7.58E-5,  2.32E-3])

for osc in oscparams:

    update = "false"
    for mode in modes:

        command = 'root2k -l -b -q \'$T2KOSCANA/macros/physin/test_sk_ereco_spectrum_twk_sys.C('+str(mode)+','+update+',"3f_osc",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+','+str(emax)+',"'+conffile+'")\''
        print command
        if os.system(command) == 1:
            sys.exit(1)
        if update == "false":
            update = "true"

"""
new TCanvas();
categories_p->Draw();
categories_m->Draw("same");
categories_p->GetYaxis()->SetRangeUser(-100,100);
categories_p->SetTitle("");
categories_p->GetXaxis()->SetLabelSize(0.03); 
categories_m->SetLineColor(kBlue);
TLine line_c(0,0,categories_p->GetNbinsX(),0);
line_c->Draw();
new TCanvas();
summary_p->Draw();
summary_m->Draw("same");
summary_p->GetYaxis()->SetRangeUser(-10,10);
summary_p->SetTitle("");
summary_p->GetXaxis()->SetLabelSize(0.03); 
summary_m->SetLineColor(kBlue);
TLine line_s(0,0,summary_p->GetNbinsX(),0);
line_s->Draw();

"""
