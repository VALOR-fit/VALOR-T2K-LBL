#!/usr/bin/python
import ROOT as R
import sys
import glob
import os

# Alternative to 'hadd' because it's segfaulting for me.
# Only works on TTrees


# Need argument parser here for experiment - for now check if dune defined
if "VALOR_DUNE" in os.environ:
  experiment = "DUNE"
else:
  experiment = "T2K"
    
try:
  ROOT.valor
except:
  if   experiment=="T2K" : R.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C")
  elif experiment=="DUNE": R.gROOT.ProcessLine(".x $VALOR_DUNE/root/valor_dune.C")
  else:                    R.gROOT.ProcessLine(".x $VALOR_FMWK/root/valor_interpreter.C");

output_file = sys.argv[1]
input_files = sys.argv[2:]

first_file = glob.glob(input_files[0])[0]

ff = R.TFile(first_file)

tree_names = [k.GetName() for k in ff.GetListOfKeys() if k.ReadObj().ClassName()=="TTree"]

ff.Close()

of = R.TFile(output_file,"CREATE")
if not of.IsOpen():
  print "Failed to open output file ",of," - does it already exist? This script won't overwrite."
  exit()
  
of.cd()

for tree_name in tree_names:
  chain = R.TChain(tree_name)
  for f in input_files:
    chain.Add(f)
  
  output_tree = chain.CloneTree()
  of.cd()
  output_tree.Write()

of.Close()
