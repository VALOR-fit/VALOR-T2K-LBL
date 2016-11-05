#!/usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging
import glob

"""
This script will run the macro test_sk_ereco_spectrum_twk_sys_throw_nuebar.C for all systematic possiblities, producing
a ROOT file with the outputs. Things you can change are:
 - postfit_config - configuration file for BANFF postfit
 - prefit_config - configuration file for BANFF prefit
 - oscparam - the oscillation parameters you want to calculate shifts relative to (defaults to global best fits - you can specify more than 1 set of parameters to run over)
            - Give each of the oscparam lists a tag and they will be compared side by side
            - Also specify which oscillation parameters are to be randomized by a string eg ss23ss13dcpss12dmsq23dmsq12 or ss23ss13dcpss12 (Macro uses reactor constraint by default)
 - save_toy_expts - save the toy experiment ereco spectra? Warning: will produce big files! (20000 toy spectra per syst!)
 - make_summary - should we write a summary table for publication
 
Run this by /path/to/test_sk_ereco_spectrum_twk_sys_throw_runAll.py [--all] (you may need to chmod +x)
[if the --all option is set, will run systematics individually also, on top of the normal 5 groups]
It will then create the .root in the current working directory
"""

#
# NEW CONFIGURATION FILES (4/6/16 BANFF RELEASE W/ FIT BIAS AND ADDITIONAL UNCERTAINTY ON RUN 1-7 FLUX)
#
#postfit_config  = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17.xml"
#prefit_config  = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FluxTunRun17_prefit.xml"

#
# NEW CONFIGURATION FILES (14/3/16 BANFF RELEASE W/ AND W/O FIT BIAS)
#
#postfit_config  = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias.xml"
#postfit_config  = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_FSIOnly.xml"
#postfit_config  = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_SKDetOnly.xml"
#postfit_config  = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias.xml"
#postfit_config  = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias_ExpSummer2016.xml"
#prefit_config  = "analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_fitbias_prefit.xml"

prior_file = ""

table_fn = "nsk_systematic_table_summary_Asimov_A.tex"

#osc_tag_file = "dcp_0.000" # Asimov B
osc_tag_file = "dcp_-99.000" # Asimov A

unosc_tag_file = "dcp_0.000"

class T2KAnaError(Exception): pass

optionParser = optparse.OptionParser(usage ="""\
usage: %prog [options]
    Run test_sk_ereco_spectrum_twk_sys_throw_joint_2015_2d.C many times""")
#names of things
optionParser.add_option("-a","--all", action="store_true", dest="do_all", default=False, help="If set, will run indiviudal systematics, as well as the categories")
optionParser.add_option("-i","--important", action="store_true", dest="do_imp", default=False, help="If set, will run the spline/osc param indiviudal systematics, as well as the categories")
optionParser.add_option("-n","--nexp", action="store", dest="nexp", type="int", default=-1, help="How many experiments to generate at each point. Required")
optionParser.add_option("--no-table", action="store_true", dest="notable", default=False, help="Don't make a table summarising the results for a tech note")
optionParser.add_option("-s","--sksample", action="store", dest="sksample", type="string", default="1Rmu", help="Choose the SK sample: 1Rmu, 1Re, RHC1Rmu, RHC1Re")
optionParser.add_option("-e","--hierarchy", action="store", dest="hierarchy", type="string", default="normal", help="Mass Hierarchy: normal, inverted")
optionParser.add_option("-m","--mixparamchoice", action="store", dest="mixparamchoice", type="string", default="sinsqtheta", help="Mixing param choice: sinsqtheta")
optionParser.add_option("-f","--fitparamlist", action="store", dest="fitparamlist", type="string", default="nuosc3f_12fix+nuisance", help="Mixing param choice: sinsqtheta")
optionParser.add_option("--hyp","--hypothesis", action="store", dest="hypothesis", type="string", default="3f_osc", help="Physics hypothesis")
optionParser.add_option("-t","--storethrows", action="store_true", dest="storethrows", default=False, help="Store the thrown parameter values")
optionParser.add_option("--do-not-run", action="store_true", dest="do_not_run", default=False, help="If you've already got the files and just want to use them rather than running again.")
optionParser.add_option("--mode",action="store",dest="mode",default=None,help="If you only want to run a single mode, specify its number here")
optionParser.add_option("--skip-prefit" , action="store_true", dest="skip_prefit" , default=False, help="If don't want to run prefit jobs")
optionParser.add_option("--skip-postfit", action="store_true", dest="skip_postfit", default=False, help="If don't want to run postfit jobs")
optionParser.add_option("--skip-osc"  , action="store_true", dest="skip_osc"  , default=False, help="If don't want to run oscillated jobs")
optionParser.add_option("--skip-unosc", action="store_true", dest="skip_unosc", default=False, help="If don't want to run unoscillated jobs")
optionParser.add_option("--randomise-osc", action="store", dest="random_osc", default="ss12-dmsq12-ss13-", help="Which params I want randomized in random osc toys, must have a - at the end!")
optionParser.add_option("--clear", action="store_true", dest="clear", default=False, help="Get rid of output of other jobs first.")
optionParser.add_option("--config-postfit", action="store", dest="postfit_config", default="analysis_configuration_data_theta_1bin_NIWG2015c_1p1h_fitbias_Run17c_final.xml")
optionParser.add_option("--config-prefit", action="store", dest="prefit_config", default="analysis_configuration_data_theta_1bin_NIWG2015c_1p1h_fitbias_Run17c_final_prefit.xml")

optionParser.add_option("-p","--prob-nh", action="store", dest="prob_nh", default=0.5, help="Prior for Normal Hierarchy")
optionParser.add_option("-d","--do-marg", action="store", dest="domarg", default=False, help="Flag for marginalization")
#optionParser.add_option("-x","--mixparamoctant", action="store", dest="mixparamoctant", type="string", default="two", help="Mixing param octant: two")
#optionParser.add_option("-z","--rethrow-unphysical-toys", action="store", dest="rethrow_toys", default=False, help="Flag to rethrow unphysical toys")
optionParser.add_option("-u","--dnovern-new-def", action="store", dest="dnovern", default=True, help="Use new definition for dn/n")

#optionParser.add_option("-z","--do-totsyst-squaresum", action="store", dest="dototsyssum", default=True, help="")
#optionParser.add_option("-p","--do-plot-osc", action="store", dest="plotosc", default=True, help="")

(opt, args) = optionParser.parse_args()

make_summary = not opt.notable
postfit_config = opt.postfit_config
prefit_config = opt.prefit_config
prob_nh = opt.prob_nh
if prob_nh < 0: 
    sys.exit("Prob(NH) must be positive")
if prob_nh > 1: 
    sys.exit("Prob(NH) must be below 1")
dnovern = "true" if opt.dnovern else "false"
mixparamoctant = "two"
#mixparamoctant = opt.mixparamoctant
rethrow_toys = "true"
#rethrow_toys = opt.rethrow_toys
domarg = "true" if opt.domarg else "false"
#domarg = opt.domarg
dmsqscheme = "auto"

nexp = opt.nexp
if nexp <= 0:
    sys.exit("specify number of events")
sksample = opt.sksample
hierarchy = opt.hierarchy
mixparamchoice = opt.mixparamchoice
hypothesis = opt.hypothesis
fitparamlist = opt.fitparamlist
do_not_run_jobs = opt.do_not_run
store_throws = "true" if opt.storethrows else "false"
random_osc = opt.random_osc
if args:
  print "Unrecognised arguments: ",args
  exit(-1)

save_toy_expts = "false"

#   0 tweak all parameters (name_modes="all")
#   1 tweak all BANFF (name_modes="BANFF")
#   5 tweak all BANFF oxygen (name_modes="") --> not any longer 
#   6 tweak all BANFF non-oxygen (name_modes="") --> not any longer
#   7 tweak BANFF flux only (name_modes="BANFF_flux")
#   8 tweak BANFF xsec only (name_modes="BANFF_xsec")
#   9 tweak SKdet + FSI (name_modes="SKDet_FSI")
#  11 oscillation parameters (name_modes="Osc_param")
#  12 tweak nuenumu (nu e nubar)
#  13 tweak SKDetFSI+BANFF (all except oscillations)
# 650 tweak mec fhc+rhc (name_modes="banff50_fhcprhc")
# 666 tweak mecbar fhc+rhc (name_modes="banff66_fhcprhc")
# 662 tweak NC other (name_modes="banff62_fhcprhc")
# 664 tweak NC 1gamma (name_modes="banff64_fhcprhc")

#
# Set the systematic running mode 
# Keep the order you want to have in the table
# Keep same order between modes and name_modes
#

name_modes = []
latex_name_modes = []
modes = []

# modes = 9
modes.append(9)
name_modes.append("SKDet_FSI") 
latex_name_modes.append("SKDet+FSI+SI")
# modes = 7
modes.append(7)
name_modes.append("BANFF_flux") 
latex_name_modes.append("Flux")
# modes = 8
modes.append(8)
name_modes.append("BANFF_xsec") 
latex_name_modes.append("XSec")
# modes = 650
modes.append(650)
name_modes.append("banff50_fhcprhc")
latex_name_modes.append("MEC")
# modes = 666
modes.append(666)
name_modes.append("banff66_fhcprhc")
latex_name_modes.append("MECbar")
# modes = 662
modes.append(662)
name_modes.append("banff62_fhcprhc")
latex_name_modes.append("NC other")
# modes = 664
modes.append(664)
name_modes.append("banff64_fhcprhc")
latex_name_modes.append("NC 1gamma")
# modes = 12
modes.append(12)
name_modes.append("NueNumu") 
latex_name_modes.append("XSec nue/numu")
# modes = 13
#modes.append(13)
#name_modes.append("skdetfsi_banff") 
#latex_name_modes.append("Flux+BANFF+SKDetFSI+SI")
# modes = 1
modes.append(1)
name_modes.append("BANFF") 
latex_name_modes.append("BANFF Tot")
# modes = 1 (prefit)
modes.append(1)
name_modes.append("BANFFprefit") 
latex_name_modes.append("BANFF Tot (prefit)")
# modes = 11
modes.append(11)
name_modes.append("Osc_param") 
latex_name_modes.append("Oscillations")
# modes = 0
modes.append(0)
name_modes.append("all") 
latex_name_modes.append("All")
# modes = 0 (prefit)
modes.append(0)
name_modes.append("allprefit") 
latex_name_modes.append("All (prefit)")



if opt.mode is not None:
  modes = [int(opt.mode)]

oscparams = [] #fill with oscillation parameters, in order mixing12, theta13, theta23, deltacp, dmsq12, dmsq23, betabar hypothesis
oscparams_tag = []

if not opt.skip_osc  :
  oscparams.append([0.304, 0.02172183825, 0.528, -1.601, 0.0000753, 0.002509,  -99]) # use global best fits!
  oscparams_tag.append("Toy_A")

#if not opt.skip_unosc: oscparams.append([  0,   0,   0,   0,   -99,   -99,  -99]) #no oscillations!

# if not opt.skip_osc  :
#   oscparams.append([0.304, 0.0217, 0.45,  0, 0.0000753, 0.002509,  -99]) # use global best fits!
#   oscparams_tag.append("Toy_B")

#if not opt.skip_osc  : oscparams.append([-99, -99, -99, -99, -99, -99,  -99]) # use global best fits!
#if not opt.skip_osc  : oscparams.append([-999, -999, -999, -99, -999, -999, -99]) # fix dcp
#if not opt.skip_osc  : oscparams.append([-999, -999, -99, -99, -999, -99, -99]) # fix theta23, dmsq23 and dcp

update = "true"

if not do_not_run_jobs:
  for osc_index, osc in enumerate(oscparams,start=0):

      if opt.clear: update = "false"

      for mode in modes:

          command = 'root2k -l -b -q \'$T2KOSCANA/macros/spectra/test_sk_ereco_spectrum_twk_sys_throw_joint_2015_2d.C('+str(mode)+','+update+',"'+sksample+'","'+hypothesis+'",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+','+str(osc[6])+','+save_toy_expts+','+store_throws+',"'+str(postfit_config)+'",'+str(nexp)+',"'+dmsqscheme+'","'+mixparamchoice+'","'+hierarchy+'","'+mixparamoctant+'",'+rethrow_toys+','+dnovern+','+domarg+','+str(prob_nh)+',"'+prior_file+'","'+oscparams_tag[osc_index]+'","'+random_osc+'")\''


          print command
          if not opt.skip_postfit:
            if os.system(command) == 1:
                sys.exit(1)
            if update == "false": 
                update = "true"
        
          
          #repeat for the banff prefit (its a different spectrum, so erros have different effects)
          if not opt.skip_prefit:
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/spectra/test_sk_ereco_spectrum_twk_sys_throw_joint_2015_2d.C('+str(mode)+','+update+',"'+sksample+'","'+hypothesis+'",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+','+str(osc[6])+','+save_toy_expts+','+store_throws+',"'+str(prefit_config)+'",'+str(nexp)+',"'+dmsqscheme+'","'+mixparamchoice+'","'+hierarchy+'","'+mixparamoctant+'",'+rethrow_toys+','+dnovern+','+domarg+','+str(prob_nh)+',"'+prior_file+'","'+oscparams_tag[osc_index]+'","'+random_osc+'")\''

            print command
            if os.system(command) == 1:
                sys.exit(1)


if make_summary:
  # Synthesise nsk summary table
  # Eg table 6 in TN-155
  # Uses outputs from previous macro runs, but calculates 'all' values using quadrature.
  files = glob.glob("*_*txt")
    
  col_titles = []
  col_files = []
  
  if "RHC" not in sksample :
    p = [x for x in files if sksample in x and "RHC" not in x ]
  else:
    p = [x for x in files if sksample in x ]

  sample_info = p[0]

  for osctag in oscparams_tag:
    l = [x for x in p if osctag in x ]
    #l = [x for x in files if osc_tag_file in x ]
    if not len(l)==1:
      raise T2KAnaError("Found "+`len(l)`+ "txt files for "+osctag+", expected 1")
    col_files.append(l[0])
    col_titles.append(osctag)


  #col_titles = ["(no oscillations)","(oscillations with typical parameter values)"]
  #col_files = [no_osc_info,osc_info]


  out_cols = []
  for title,infile in zip(col_titles,col_files):
    f = open(infile)
    input_rows = dict()
    for l in f:
      try:
        row_id,val = l.split("& &")
        row_id = row_id.strip()
        val = val.strip().partition(" ")[0]
        input_rows[row_id] = val

      except ValueError:
        pass
    f.close()
    d = dict()
    out_cols.append(d)


    for name in name_modes:
        d[name] = float(input_rows[name])
            
    
    # Add XSec Corr Tot (TO CHECK)
    if "BANFF" in input_rows and "banff62_fhcprhc" in input_rows and "banff64_fhcprhc" in input_rows and "NueNumu" in input_rows:
      name_curr = "xsec_corr"
      name_modes.append(name_curr)
      latex_name_modes.append("XSec Corr Tot")
      d[name_curr] = math.sqrt( float( float(input_rows["BANFF"])**2 - float(input_rows["banff62_fhcprhc"])**2 - float(input_rows["banff64_fhcprhc"])**2 - float(input_rows["NueNumu"])**2 ) )


    # Add Flux+XSec+SKDetFSI
    if "BANFF" in input_rows and "SKDet_FSI" in input_rows:
      name_curr = "flux_xsec_skdetfsi_calc"
      name_modes.append(name_curr)
      latex_name_modes.append("Flux+XSec+SKDet+FSI+SI (Calculated)")
      d[name_curr] = math.sqrt( float( float(input_rows["BANFF"])**2 + float(input_rows["SKDet_FSI"])**2 ) )

    

  f = open(table_fn,"a") #append new table to file.
  f.write("""\\begin{table}[htb]
      \\centering
      \\begin{tabular}{| c | c | c |}
      \\hline\n""")

  # Takes a list of column contents (as strings), makes it into a latex table row, roughly fixed width
  def make_row(list_of_columns,indent=2,width=30):
    fields = [x.ljust(width) for x in list_of_columns]
    return (" "*indent) + (" & ".join(fields))

  # Formats a decimal into a percentage, to 2dp, for latex
  def percent_2dp(v):
    return "{0:.2f}".format(100*v)+"\\%"
    
  f.write(make_row(["{\\bf Source of uncertainty}"]+(["\\multicolumn{1}{c|}{{\\bf ${\\delta}N_{SK}$/$N_{SK}$}}"]*len(col_titles))))
  f.write("\\\\\n")
  f.write(make_row([""]+["\multicolumn{1}{c|}{"+t+"}" for t in col_titles]))
  f.write("\\\\\n\\hline\n")

  #
  #for latex_name,name in zip(latex_name_modes,name_modes):
  #    print "\n",latex_name," - ",name,"\n"
  #    f.write(make_row([latex_name]+[percent_2dp(d[name]) for d in out_cols]))
  #    f.write("\\\\\n")
  #
 
  namerow = "SKDet_FSI"
  if namerow in name_modes:
      f.write("\\hline")
      f.write(make_row(["SKDet+FSI+SI"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      f.write("\\hline\n")

  namerow = "BANFF_flux"
  if namerow in name_modes:
      f.write(make_row(["Flux"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      f.write("\\hline\n")

  namerow = "xsec_corr"
  if namerow in name_modes:
      f.write(make_row(["XSec Corr Tot"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      #f.write("\\hline\n")

  namerow = "banff50_fhcprhc"
  if namerow in name_modes:
      f.write(make_row(["MEC (corr)"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      #f.write("\\hline\n")

  namerow = "banff66_fhcprhc"
  if namerow in name_modes:
      f.write(make_row(["MECbar (corr)"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      #f.write("\\hline\n")

  namerow = "banff62_fhcprhc"
  if namerow in name_modes:
      f.write(make_row(["NC other (uncorr)"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      #f.write("\\hline\n")

  namerow = "banff64_fhcprhc"
  if namerow in name_modes:
      f.write(make_row(["NC 1gamma (uncorr)"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      #f.write("\\hline\n")

  namerow = "NueNumu"
  if namerow in name_modes:
      f.write(make_row(["XSec nue/numu"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      #f.write("\\hline\n")

  namerow = "BANFF_xsec"
  if namerow in name_modes:
      f.write(make_row(["XSec Tot"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      f.write("\\hline\n")

  namerow = "BANFF"
  if namerow in name_modes:
      f.write(make_row(["Flux+XSec"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      f.write("\\hline\n")

  namerow = "skdetfsi_banff"
  if namerow in name_modes:
      f.write(make_row(["Flux+XSec+SKDet+FSI+SI"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      f.write("\\hline\n")

  namerow = "flux_xsec_skdetfsi"
  if namerow in name_modes:
      f.write(make_row(["Flux+XSec+SKDet+FSI+SI"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      f.write("\\hline\n")
      
  namerow = "Osc_param"
  if namerow in name_modes:
      f.write(make_row(["Oscillations"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      f.write("\\hline\n")

  namerow = "all"
  if namerow in name_modes:
      f.write(make_row(["All"]+[percent_2dp(d[namerow]) for d in out_cols]))
      f.write("\\\\\n")
      f.write("\\hline\n")

  f.write("\\hline\n")
  
  f.write("""\\end{tabular}
  %caption goes here
  \\caption{""")
  f.write(sksample)
  f.write("""}
  \\label{tab:nsk_systematic_table_summary}
  \\end{table}\n""")

  f.close()
  
  print "Printed LaTeX table to ",table_fn

