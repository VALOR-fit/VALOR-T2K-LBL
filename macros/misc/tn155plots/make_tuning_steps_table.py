#!/usr/bin/env python

# Script to make table of nominal/flux-tuned/banff-tuned spectra
#
# Requires you to of already created the output of
# test_sk_ereco_spectrum_nominal.C 
# (since this needs to be run & recompiled a few times)

import ROOT
import math
import os
import optparse

#program = "root2k"
#prog_args = ["-q","-b"]
ana_dir = os.environ["T2KOSCANA"]
try:
  ROOT.t2kana
except:
  ROOT.gROOT.ProcessLine(".x "+ana_dir+"/root/loadincs.C");    
  ROOT.gROOT.ProcessLine(".x "+ana_dir+"/root/loadlibs.C");
ROOT.gROOT.SetBatch(1)
conv = ROOT.t2kana.conv
util = ROOT.t2kana.util

optionParser = optparse.OptionParser(usage ="""\
usage: %prog [options]""")
optionParser.add_option("--osc", action="store_true", dest="osc", default=False, help="If set, will do oscillated table")
(opt, args) = optionParser.parse_args()
osc = opt.osc
print osc

if osc:
    with open('banff_tuned/sk_ereco_spectrum-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2.txt', 'r') as f_banfftuned:
        a_banfftuned = f_banfftuned.readlines()
    f_banfftuned.close()
    with open('flux_tuned/sk_ereco_spectrum-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2.txt', 'r') as f_fluxtuned:
        a_fluxtuned = f_fluxtuned.readlines()
    f_fluxtuned.close()
    with open('nominal/sk_ereco_spectrum-sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2.txt', 'r') as f_nominal:
        a_nominal = f_nominal.readlines()
    f_nominal.close()
    f = open('nsk_nominal_osc.tex','w')
else:
    with open('banff_tuned/sk_ereco_spectrum-sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00000000,23_0.000000]eV2.txt', 'r') as f_banfftuned:
        a_banfftuned = f_banfftuned.readlines()
    f_banfftuned.close()
    with open('flux_tuned/sk_ereco_spectrum-sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00000000,23_0.000000]eV2.txt', 'r') as f_fluxtuned:
        a_fluxtuned = f_fluxtuned.readlines()
    f_fluxtuned.close()
    with open('nominal/sk_ereco_spectrum-sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00000000,23_0.000000]eV2.txt', 'r') as f_nominal:
        a_nominal = f_nominal.readlines()
    f_nominal.close()
    f = open('nsk_nominal_noosc.tex','w')


#loop over the lines
for i, value in enumerate(a_nominal):
    banff   = a_banfftuned[i]
    flux    = a_fluxtuned[i]
    nominal = a_nominal[i]
    #find out if it's the first (total) line
    mode = nominal.split('otal number of events: ')
    if len(mode) == 2:
        #header stuff
        f.write('\\begin{table}\n')
        f.write('\\centering\n')
        f.write('\\begin{tabular}{| c |  D{.}{.}{-1} |  D{.}{.}{-1} |  D{.}{.}{-1} | }\n')
        f.write('\\hline\n')
        f.write('         & \multicolumn{1}{c|}{$N_{SK}$}     & \multicolumn{1}{c|}{$N_{SK}$}          & \multicolumn{1}{c|}{$N_{SK}$}          \\\\ \n')
        f.write('         & \multicolumn{1}{c|}{(nominal MC)} & \multicolumn{1}{c|}{(NA61-tuned MC)}  & \multicolumn{1}{c|}{(BANFF-tuned MC)}  \\\\ \n')
        f.write('\\hline\n')
        #total line
        nominal_this = nominal.split('otal number of events: ')[1].strip()
        flux_this    = flux.split   ('otal number of events: ')[1].strip()
        banff_this   = banff.split  ('otal number of events: ')[1].strip()
        f.write('Total & '+nominal_this+' & '+flux_this+' & '+banff_this+' \\\\ \n')
        f.write('\\hline\n')
        continue
    #find out if it's an individual mode line
    mode = nominal.split('Contrib from mode: ')
    if len(mode) == 1:
        continue
    mode_this  = nominal.split('Contrib from mode: ')[1].split(':')[0].strip()
    superkmode = util.SuperkModeFromString(mode_this)
    #print util.AsStringTex.__doc__
    mode_this  = util.AsStringTex.disp('t2kana::conv::SuperkMode_t m')(superkmode)
    mode_this  = mode_this.split('(')[0].strip() #get rid of '(numu disap)'
    nominal_this = nominal.split('Contrib from mode: ')[1].split(':')[1].split('events')[0].strip()
    flux_this    = flux.split   ('Contrib from mode: ')[1].split(':')[1].split('events')[0].strip()
    banff_this   = banff.split  ('Contrib from mode: ')[1].split(':')[1].split('events')[0].strip()
    f.write(mode_this+' & '+nominal_this+' & '+flux_this+' & '+banff_this+' \\\\ \n')
#footer
f.write('\\hline\n')
f.write('\\end{tabular}\n')
f.write('\\caption{\n')
if osc:
    f.write('Calculated numbers of 1-ring $\\mu$-like events using the nominal MC templates, the NA61-tuned MC templates\n')
    f.write('and the BANFF-tuned MC templates.\n')
    f.write('The total numbers of events and the numbers of events from each mode considered in this analysis are shown.\n')
    f.write('These numbers were calculated for an exposure of 6.5700$\\times$10$^{20}$ POT\n')
    f.write('and oscillations with $\\sin^{2}\\theta_{23}$=0.5 and $|\\Delta m^{2}_{32}|$=2.4$\\times$10$^{-3}$ eV$^{2}$/c$^{4}$.\n')
    f.write('Values shown in Tab. \\ref{tab:pdg2012non23} are used for $\\sin^{2}\\theta_{13}$, $\\sin^{2}\\theta_{12}$ and $\\Delta m^{2}_{21}$ and $\\delta_{CP}$=0.\n')
    f.write('The normal mass hierarchy is assumed.\n')
    f.write('}\n')
    f.write('\\label{tab:nsk_nominal_x3_default_osc}\n')
else:
    f.write('  Calculated numbers of 1-ring $\\mu$-like events without oscillations using the nominal MC templates, the\n')
    f.write('  NA61-tuned MC templates and the BANFF-tuned MC templates.\n')
    f.write('  The total numbers of events and the numbers of events from each mode considered in this analysis are shown.\n')
    f.write('  These numbers were calculated for an exposure of 6.5700$\\times$10$^{20}$ POT.\n')
    f.write('}\n')
    f.write('\\label{tab:nsk_nominal_x3_no_osc}\n')
f.write('\\end{table}\n')

