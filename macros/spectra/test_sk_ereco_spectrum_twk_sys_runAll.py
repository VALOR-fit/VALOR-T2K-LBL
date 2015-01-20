#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging
import glob

"""
This script will run the macro test_sk_ereco_spectrum_twk_sys.C for all systematic possiblities, producing
a ROOT file with the outputs. Things you can change are:
 - postfit_config - configuration file for BANFF postfit
 - prefit_config - configuration file for BANFF prefit
 - oscparam - the oscillation parameters you want to calculate shifts relative to (defaults to global best fits - you can specify more than 1 set of parameters to run over)
            - (Note that if sinsqtheta13 is identical in 2 cases, you'll need to add another distinguisher to the label definition at the end of root macro, otherwise the tex won't compile)
 - save_toy_expts - save the toy experiment ereco spectra? Warning: will produce big files! (20000 toy spectra per syst!)
 - make_summary - should we write a summary table for publication
 
Run this by /path/to/test_sk_ereco_spectrum_twk_sys_runAll.py [--all] (you may need to chmod +x)
[if the --all option is set, will run systematics individually also, on top of the normal 5 groups]
It will then create the .root in the current working directory
"""
    
postfit_config  = "analysis_configuration.xml"

prefit_config  = "analysis_configuration.xml"


class T2KAnaError(Exception): pass

optionParser = optparse.OptionParser(usage ="""\
usage: %prog [options]
    Run test_sk_ereco_spectrum_twk_sys.C many times""")
#names of things
optionParser.add_option("-a","--all", action="store_true", dest="do_all", default=False, help="If set, will run indiviudal systematics, as well as the categories")
optionParser.add_option("-i","--important", action="store_true", dest="do_imp", default=False, help="If set, will run the spline/osc param indiviudal systematics, as well as the categories")
optionParser.add_option("--no-table", action="store_true", dest="notable", default=False, help="Don't make a table summarising the results for a tech note")
optionParser.add_option("-s","--sksample", action="store", dest="sksample", type="string", default="1Rmu", help="Choose the SK sample: 1Rmu, 1Re, RHC1Rmu, RHC1Re")
optionParser.add_option("-e","--hierarchy", action="store", dest="hierarchy", type="string", default="normal", help="Mass Hierarchy: normal, inverted")
optionParser.add_option("-m","--mixparamchoice", action="store", dest="mixparamchoice", type="string", default="sinsqtheta", help="Mixing param choice: sinsqtheta")
optionParser.add_option("-f","--fitparamlist", action="store", dest="fitparamlist", type="string", default="nuosc3f_12fix+nuisance", help="Mixing param choice: sinsqtheta")

(opt, args) = optionParser.parse_args()

make_summary = not opt.notable

sksample = opt.sksample
hierarchy = opt.hierarchy
mixparamchoice = opt.mixparamchoice
fitparamlist = opt.fitparamlist

save_toy_expts = "false"

skip_prefit = False

modes1 = [0,1,2,9,11] #combinations

#modes2 = range(400,422) #banff fhc (niwg2013 - numudisap)
#modes2 = range(500,422) #banff rhc (niwg2013 - numudisap)
#modes2 = range(500,527) #banff fhc (niwg2015 - numudisap)
modes2 = range(500,537) #banff rhc (niwg2015 - numudisap)
#~ modes2 = range(525,537) #banff rhc (niwg2015 - numudisap)

#modes3 = range(200,210) #unconstranied xsec (niwg2013)
modes3 = range(211,211) #unconstranied xsec (niwg2015)

#modes4 = range(800,805) #sk detector (skdetfsi 2013 - numudisap - fhc)
#modes4 = range(806,817) #sk detector (skdetfsi 2013 - nue - fhc)
#modes4 = range(800,817) #sk detector (skdetfsi 2013 - joint - fhc)
modes4 = range(900,905) #sk detector (skdetfsi 2013 - numudisap - rhc)
#modes4 = range(906,917) #sk detector (skdetfsi 2013 - nue - rhc)
#modes4 = range(900,917) #sk detector (skdetfsi 2013 - joint - rhc)

modes6 = range(700,703) #osc param

modes7 = [7,8] #banff flux / banff xsec individually

modes = modes1
if opt.do_all: #combinations + everything
    modes = modes + modes2 + modes3 + modes4 + modes6 + modes7
elif opt.do_imp: #combinations + splines + oscparams
    modes = modes + modes6 + [200,201,202,203,204,215]
    
    
#DEBUG!
modes = modes2 + modes4 + modes6
    

oscparams = [] #fill with oscillation parameters, in order mixing12, theta13, theta23, deltacp, dmsq12, dmsq23hypothesis

oscparams.append([  0,   0,   0,   0,   -99,   -99]) #no oscillations!
oscparams.append([-99, -99, -99, -99, -99, -99]) #use global best fits!
#oscparams.append([0.86, 0.1,  1,  0,  7.58E-5,  2.32E-3])
#oscparams.append([0.86, 0,    1,  0,  7.58E-5,  2.32E-3])

for osc in oscparams:

    update = "false"

    for mode in modes:

        #changing osc params doesn't affect anything if there are no oscillations
        if mode == 11 and osc == [0,   0,   0,   0,   0,   0]:
            continue
        if mode == 11 and osc == [0,   0,   0,   0,   -99,  -99]:
            continue
        
        command = 'root2k -l -b -q \'$T2KOSCANA/macros/spectra/test_sk_ereco_spectrum_twk_sys.C('+str(mode)+','+update+',"'+sksample+'","3f_osc",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+','+save_toy_expts+',"'+postfit_config+'","'+mixparamchoice+'","'+hierarchy+'")\''
        print command
        if os.system(command) == 1:
            sys.exit(1)
        if update == "false":
            update = "true"
            
        #repeat for the banff prefit (its a different spectrum, so erros have different effects)
        
        if not skip_prefit:
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/spectra/test_sk_ereco_spectrum_twk_sys.C('+str(mode)+','+update+',"'+sksample+'","3f_osc",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+','+save_toy_expts+',"'+prefit_config+'","'+mixparamchoice+'","'+hierarchy+'")\''
          print command
          if os.system(command) == 1:
              sys.exit(1)


if make_summary:
  # Synthesise nsk summary table
  # Eg table 6 in TN-155
  # Uses outputs from previous macro runs, but calculates 'all' values using quadrature.
  files = glob.glob("sk_ereco_spectrum_*-twk_sys_*txt")
  l = [x for x in files if "12_0.000,13_0.000,23_0.000" in x]
  if not len(l)==1:
    raise T2KAnaError("Found "+`len(l)`+ "txt files for no oscillation, expected 1")
  no_osc_info = l[0]
  l = [x for x in files if "12_0.000,13_0.000,23_0.000" not in x]
  if not len(l)==1:
    raise T2KAnaError("Found "+`len(l)`+ "txt files for oscillation, expected 1")
  osc_info = l[0]

  col_titles = ["(no oscillations)","(oscillations with typical parameter values)"]
  col_files = [no_osc_info,osc_info]

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

    d["skdet"]      = float(input_rows["SKDet_FSI"])
    d["banff_pre"]  = float(input_rows["BANFFprefit"])
    d["banff_post"] = float(input_rows["BANFF"])
    d["xsec"]       = float(input_rows["unconstrained_xsec"])
    
    all_postfit = ["SKDet_FSI","BANFF","unconstrained_xsec"]
    all_with_banff_prefit = ["SKDet_FSI","BANFFprefit","unconstrained_xsec"]

    if "no oscillations" not in title:
      all_postfit.append("Osc_param")
      all_with_banff_prefit.append("Osc_param")
      
    all_prefit = [x+"prefit" for x in all_postfit]
    
    # All from quadrature
    #~ d["total_prefit"]  = math.sqrt(sum([float(input_rows[x])**2 for x in all_prefit ])) # all using the 'prefit' figures
    d["total_prefit"]  = math.sqrt(sum([float(input_rows[x])**2 for x in all_with_banff_prefit ])) # Only BANFF prefit figures used
    d["total_postfit"] = math.sqrt(sum([float(input_rows[x])**2 for x in all_postfit]))
    
    try:
      d["osc"] = float(input_rows["Osc_param"])
    except KeyError:
      d["osc"] = 0.
      
  f = open("nsk_systematic_table_summary.txt","w")
  f.write("""\\begin{table}[htb]
\\centering
\\begin{tabular}{| c | D{.}{.}{-1} | D{.}{.}{-1} |}
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

  f.write(make_row(["SuperK detector + FSI + SI"]+[percent_2dp(d["skdet"]) for d in out_cols]))
  f.write("\\\\\n")

  f.write(make_row(["BANFF (prefit)"]+[percent_2dp(d["banff_pre"]) for d in out_cols]))
  f.write("\\\\\n")

  f.write(make_row(["BANFF (postfit)"]+[percent_2dp(d["banff_post"]) for d in out_cols]))
  f.write("\\\\\n")

  f.write(make_row(["Uncorrelated XSec"]+[percent_2dp(d["xsec"]) for d in out_cols]))
  f.write("\\\\\n")

  f.write(make_row(["\emph{Non-23-sector} oscillation parameters"]+[percent_2dp(d["osc"]) for d in out_cols]))
  f.write("\\\\\n\\hline\n")

  f.write(make_row(["Total (BANFF prefit)"]+[percent_2dp(d["total_prefit"]) for d in out_cols]))
  f.write("\\\\\n")

  f.write(make_row(["Total (BANFF postfit)"]+[percent_2dp(d["total_postfit"]) for d in out_cols]))
  f.write("\\\\\n\\hline\n")
  
  f.write("""\\end{tabular}
%caption goes here
\\label{tab:nsk_systematic_table_summary}
\\end{table}\n""")


  f.close()

