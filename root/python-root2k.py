print "\n============================="
print   "Loading python with T2KOSCANA"
print   "=============================\n"

import ROOT
import sys
import os
import glob
import math
try:
  import numpy as np
except:
  print "Warning: Numpy not loaded"

ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");

ROOT.gROOT.ProcessLine("using namespace valor");
ROOT.gROOT.ProcessLine("using namespace util");
ROOT.gROOT.ProcessLine("using namespace conv");

rootargs = [x for x in sys.argv if x.endswith(".root")]
i=0
for f in rootargs:
  varname = "_file"+`i`
  print "Attaching file ",f," as ",varname
  globals()[varname] = ROOT.TFile(f)
  i+=1

R = ROOT
#~ t2kana = ROOT.t2kana
#~ conv = t2kana.conv
#~ util = t2kana.util
#~ Opt  = t2kana.Opt.Instance()
#~ Rnd  = t2kana.Rnd.Instance()

valor = ROOT.valor
conv = valor.conv
util = valor.util
opt    = valor.Opt.Instance()
env    = valor.Env.Instance()
t2kopt = valor.T2KOpt.Instance()
rnd    = valor.Rnd.Instance()
arridx = util.ArrIdx.Instance()
physlim = util.T2KPhysicalParameterLimits.Instance()
nr = util.NameRegistry.Instance()

env.SetDie(False) # We're interactive, so don't die unless you have to.

# Make float formatting strings
for n in xrange(1,9):
  globals()["to"+`n`+"dp"] = "{0:."+`n`+"f}"
  globals()["to"+`n`+"sf"] = "{0:."+`n-1`+"e}"
  
# Quick and dirty function to get every object in a file
def fileburst(f):
  for k in f.GetListOfKeys():
    globals()[k.GetName()] = k.ReadObj()
  return
