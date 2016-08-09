#!/usr/bin/python

import ROOT
import math
import array
import sys
import os

ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()
arridx= util.ArrIdx.Instance()

samples   = ["1Rmu","1Re","RHC1Rmu","RHC1Re","MultiRe"]

col_width = 30
banff_col_widths  = [12,60,8,12]
detfsi_col_widths = [12,90,10]

line_width = 1

to3dp = "{0:.3f}"

#______________________________________________________________________________________________________________________________

def main(sample,output_file,config_file,prefit_config):
  
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+prefit_config)
  banff = skexp.GetBANFF()
  
  prefit_err = ROOT.TMatrixD(banff.RHCCovMatrix())
  del skexp
  
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file)
  banff = skexp.GetBANFF()
  
  bf_vec = ROOT.TVectorD(banff.RHCBestFitValues())
  postfit_err = ROOT.TMatrixD(banff.RHCCovMatrix())
  det_matr = ROOT.TMatrixD(skexp.GetSuperkDetectorFSI().RHCCovMatrix())
  
  banff_bestfit       = dict()
  banff_prefit_error  = dict()
  banff_postfit_error = dict()
  detfsi_error        = dict()
  
  for ip in xrange(conv.kNumOfRHCParamsInBANFFfit):
    banff_bestfit      [ip] = bf_vec[ip]
    banff_prefit_error [ip] = math.sqrt(prefit_err [ip][ip])
    banff_postfit_error[ip] = math.sqrt(postfit_err[ip][ip])
  
  for ip in xrange(conv.kNumOfRHCSKDetFSIParams):
    detfsi_error[ip] = math.sqrt(det_matr[ip][ip])
  
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

""")

  # Write BANFF table
  tex_output.write("\\begin{table}\n")
  #~ tex_output.write("\\centering\n")
  tex_output.write("\n\n\\begin{tabular} { | c ||"+("c | "*3)+" | }\n\\hline\n")
  
  cols = ["Parameter","Description","Best fit","$1\sigma$ pre/postfit"]
  entries = [x.ljust(w) for x,w in zip(cols,banff_col_widths)]
  tex_output.write("& ".join(entries) + " \\\\ \n")
  cols = ["","","","fractional error"]
  entries = [x.ljust(w) for x,w in zip(cols,banff_col_widths)]
  tex_output.write("& ".join(entries) + " \\\\ \n")
  tex_output.write("\n\\hline\n")
  
  for ip in xrange(conv.kNumOfRHCParamsInBANFFfit):
    p = getRHCBANFFparam(ip)
    cols = [str(ip),util.ParamAsStringPaperTable(p),to3dp.format(banff_bestfit[ip])]
    cols.append(to3dp.format(banff_prefit_error [ip])+" / "+to3dp.format(banff_postfit_error [ip]))
    entries = [x.ljust(w) for x,w in zip(cols,banff_col_widths)]
    tex_output.write("& ".join(entries) + " \\\\ \n")
  
  tex_output.write("\n\\hline\n")
  tex_output.write("\n\\end{tabular}")
  tex_output.write("\\end{table}\n\n")
  
  # Write SKDetFSI table
  tex_output.write("\\begin{table}\n")
  #~ tex_output.write("\\centering\n")
  tex_output.write("\n\n\\begin{tabular} { | c ||"+("c | "*2)+" | }\n\\hline\n")
  
  cols = ["Parameter","Description","$1\sigma$ error"]
  entries = [x.ljust(w) for x,w in zip(cols,detfsi_col_widths)]
  tex_output.write("& ".join(entries) + " \\\\ \n")
  tex_output.write("\n\\hline\n")
  
  for ip in xrange(conv.kNumOfRHCSKDetFSIParams):
    p = getRHCDetFSIparam(ip)
    cols = [str(ip),util.ParamAsStringPaperTable(p),to3dp.format(detfsi_error[ip])]
    entries = [x.ljust(w) for x,w in zip(cols,detfsi_col_widths)]
    tex_output.write("& ".join(entries) + " \\\\ \n")
  
  tex_output.write("\n\\hline\n")
  tex_output.write("\n\\end{tabular}")
  tex_output.write("\\end{table}\n\n")
  
  
  tex_output.write("\\end{document}\n")
  
  print "Table output written to ",output_file

def getRHCBANFFparam(ip):
  for p in xrange(conv.knppmax):
    test_ip = arridx.RHCBANFFSystIdx(p)
    if test_ip == ip:
      return p
  print "Failed to get RHC BANFF param ",ip
  exit(-1)
  
def getRHCDetFSIparam(ip):
  for p in xrange(conv.knppmax):
    test_ip = arridx.RHCSKDetFSISystIdx(p)
    if test_ip == ip:
      return p
  print "Failed to get RHC SKDetFSI param ",ip
  exit(-1)
  
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of applied systematics")
  
  parser.add_option('-f', '--file', dest='output_file', default="size_of_systematics.tex", help="Name for output TeX file.")
  
  parser.add_option('-s', '--sample' , dest='sample'      , default="1Rmu"                    , help="Sample (1Rmu, 1Re, RHC1Rmu, RHC1Re)")
  
  parser.add_option('-c', '--config'        , dest='config' , default="analysis_configuration.xml", help="Configuration file")
  parser.add_option('-p', '--config-prefit' , dest='pconfig', default="analysis_configuration_prefit.xml", help="Prefit config file")
  
  
  o,other = parser.parse_args()
  sample = util.SKSampleFromString(o.sample)
  main(sample,o.output_file, o.config, o.pconfig)
  
  if other:
    print "Warning - ignored arguments: ",other
