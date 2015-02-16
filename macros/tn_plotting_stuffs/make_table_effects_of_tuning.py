#!/usr/bin/python

import ROOT
import math
import array
import sys
import os

env_var = "T2KOSCANA"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()
arridx= util.ArrIdx.Instance()

samples   = ["1Rmu","1Re","RHC1Rmu","RHC1Re"]

col_width = 30

to3dp = "{0:.3f}"

#______________________________________________________________________________________________________________________________

def main(sample,output_file,config_file,config_file_noflux,config_file_prefit):
  
  nominal_modes      = dict()
  NA61_flux_modes    = dict()
  prefit_tuned_modes = dict()
  postfit_tuned_modes  = dict()
  
  all_dicts = [nominal_modes,NA61_flux_modes,prefit_tuned_modes,postfit_tuned_modes]
  
  #  Make no flux SKExp
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file_noflux)
  
  #~ skexp.GetBANFF().SetTweaksToNominal()
  skexp.GetBANFF().RestoreNominalMC()
  skexp.GetBANFF().RestoreRHCNominalMC()
  #~ skexp.GetBANFF().RestoreFHCpRHCNominalMC()
  skexp.GetSuperkMC().FillTemplateArrays(False)
  
  #~ fparm = util.NamedParameterList("nuosc3f")
  phyp = valor.PhysHypothesis(conv.kHypVoid)
  skexp.Generate(phyp)
  
  nominal_modes["total"] = skexp.CurrPredictedErecoVArray(sample).Integral()
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    
    nominal_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral()
  
  del skexp
  
  # Now do the NA61 tune
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file)
  #~ skexp.GetBANFF().SetTweaksToNominal()
  skexp.GetBANFF().RestoreNominalMC()
  skexp.GetBANFF().RestoreRHCNominalMC()
  #~ skexp.GetBANFF().RestoreFHCpRHCNominalMC()
  skexp.GetSuperkMC().FillTemplateArrays(False)
  skexp.Generate(phyp)
  
  NA61_flux_modes["total"] = skexp.CurrPredictedErecoVArray(sample).Integral()
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    
    NA61_flux_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral()
  
  del skexp
  
  # NIWG tune
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file_prefit)
  skexp.Generate(phyp)
  
  prefit_tuned_modes["total"] = skexp.CurrPredictedErecoVArray(sample).Integral()
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    
    prefit_tuned_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral("")
  
  del skexp
  
  # postfit tune
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file)
  skexp.Generate(phyp)
  
  postfit_tuned_modes["total"] = skexp.CurrPredictedErecoVArray(sample).Integral()
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    
    postfit_tuned_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral("")
  
  del skexp
  
  # Start writing the TeX
  tex_output = open(output_file,"w")
  
  tex_output.write("""\documentclass[10pt]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{fullpage}
\\usepackage{geometry}

\\begin{document}

{
""")
  
  tex_output.write("\\begin{table}\n")
  #~ tex_output.write("\\centering\n")
  tex_output.write("\n\n\\begin{tabular} { | c ||"+("c | "*4)+" | }\n\\hline\n")
  
  cols = ["","$N_{SK}$","$N_{SK}$","$N_{SK}$","$N_{SK}$"]
  entries = [x.ljust(col_width) for x in cols]
  tex_output.write("& ".join(entries) + " \\\\ \n")
  
  cols = ["","Nominal MC","Flux-tuned MC","BANFF prefit MC","BANFF postfit MC"]
  entries = [x.ljust(col_width) for x in cols]
  tex_output.write("& ".join(entries) + " \\\\ \\hline \n")
  
  cols = ["Total"] + [to3dp.format(d["total"]) for d in all_dicts]
  entries = [x.ljust(col_width) for x in cols]
  tex_output.write("& ".join(entries) + " \\\\ \\hline \n")
  tex_output.write("\\hline \n")
  
  print all_dicts
  
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    cols = [util.ModeAsString(m)] + [to3dp.format(d[m]) for d in all_dicts]
    entries = [x.ljust(col_width) for x in cols]
    tex_output.write("& ".join(entries) + " \\\\ \\hline \n")
  
  
  tex_output.write("\n\\end{tabular}\n")
  tex_output.write("\\end{table}\n}\n")
  tex_output.write("\\end{document}\n")
  
  print "Table output written to ",output_file

if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of applied systematics")
  
  parser.add_option('-f', '--file', dest='output_file', default="effects_of_tuning.tex", help="Name for output TeX file.")
  
  parser.add_option('-s', '--sample' , dest='sample'      , default="1Rmu"                    , help="Sample (1Rmu, 1Re, RHC1Rmu, RHC1Re)")
  
  parser.add_option('-c', '--config'        , dest='config' , default="analysis_configuration.xml", help="Configuration file")
  parser.add_option('-p', '--config-prefit' , dest='pconfig', default="analysis_configuration_prefit.xml", help="Prefit config file")
  parser.add_option('-n', '--config-no-flux', dest='fconfig', default="analysis_configuration_nofluxtune.xml", help="Configuration file without flux tune")
  
  
  o,other = parser.parse_args()
  sample = util.SKSampleFromString(o.sample)
  main(sample,o.output_file, o.config, o.fconfig, o.pconfig)
  
  if other:
    print "Warning - ignored arguments: ",other
