print "\n==========================="
print "Loading python with T2KOSCANA"
print "===========================\n"

import ROOT
import sys
import os
import glob
import math
try:
  import numpy as np
except:
  print "Warning: Numpy not loaded"

ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/loadincs.C");    
ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/loadlibs.C");

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
Opt  = valor.Opt.Instance()
Rnd  = valor.Rnd.Instance()

