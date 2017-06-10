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

#______________________________________________________________________________________________________________________________

def banff_table_header(tex):
  tex.write("\\begin{tabular} { | " + ("c | " * 5) + " }\n\\hline \n")
  tex.write("Index & Parameter & Description & Best fit & $1\sigma$ pre/postfit\\\\ \n")
  tex.write(" & &  &  & fractional error\n\\hline \n")

def sk_table_header(tex):
  tex.write("\\begin{tabular} { | " + ("c | " * 4) + " }\n\\hline \n")
  tex.write("Index & Parameter & Description & $1\sigma$ fractional error\\\\ \n")

def table_footer(tex):
  tex.write("\n\\hline\n")
  tex.write("\n\\end{tabular}")

import numpy as np
def main(sample,output_file,prefit_config, postfit_config):

  skexp_pre = valor.SuperkExpectation()
  skexp_pre.Read(os.environ[env_var] + "/config/" + prefit_config)

  skexp_post = valor.SuperkExpectation()
  skexp_post.Read(os.environ[env_var] + "/config/" + postfit_config)
  
  params = valor.util.NamedParameterList("nuosc3f+nuisance") 
  N = params.GetNumberParameters()

  xsec_tex = open(output_file + "_xsec.tex","w")
  flux_tex = open(output_file + "_flux.tex","w")
  sk_tex = open(output_file + "_sk.tex","w")
  banff_table_header(xsec_tex)
  banff_table_header(flux_tex)
  sk_table_header(sk_tex)
  for ip in range(0, N):
    p = params.WhichParam(ip)
    name = valor.util.ParamAsStringPaperTable(p)
    if not params.IsNuisance(p):
        continue
    if valor.util.IsSKDetFSI(p):
        idx = skexp_pre.GetSuperkDetectorFSI().Param2Idx(p)
        err = np.sqrt(skexp_pre.GetSuperkDetectorFSI().CovMatrix()(idx,idx))
        sk_tex.write("%d & %s & %s & %.3f\\\\ \n" % (idx, name, "", err))
    elif valor.util.IsBANFF(p):
        idx = skexp_pre.GetBANFF().Param2Idx(p)
        bestfit = skexp_post.GetBANFF().InitialValues()(idx)
        err_pre = np.sqrt(skexp_pre.GetBANFF().CovMatrix()(idx,idx))
        err_post = np.sqrt(skexp_post.GetBANFF().CovMatrix()(idx,idx))
        if valor.util.IsXSec(p):
            xsec_tex.write("%d & %s & %s & %.3f & %.3f / %.3f\\\\ \n" % (idx, name, "", bestfit, err_pre, err_post))
        else:
            flux_tex.write("%d & %s & %s & %.3f & %.3f / %.3f\\\\ \n" % (idx, name, "", bestfit, err_pre, err_post))
  table_footer(xsec_tex)
  table_footer(flux_tex)
  table_footer(sk_tex)
  xsec_tex.close()
  flux_tex.close()
  sk_tex.close()
  
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of applied systematics")
  
  parser.add_option('-f', '--file', dest='output_file', default="size_of_systematics.tex", help="Name for output TeX file.")
  
  parser.add_option('-s', '--sample' , dest='sample'      , default="1Rmu"                    , help="Sample (1Rmu, 1Re, RHC1Rmu, RHC1Re)")
  
  parser.add_option('-c', '--config'        , dest='config' , default="analysis_configuration.xml", help="Configuration file")
  parser.add_option('-p', '--config-prefit' , dest='pconfig', default="analysis_configuration_prefit.xml", help="Prefit config file")
  
  
  o,other = parser.parse_args()
  sample = util.SKSampleFromString(o.sample)
  main(sample,o.output_file, o.pconfig, o.config)
  
  if other:
    print "Warning - ignored arguments: ",other
