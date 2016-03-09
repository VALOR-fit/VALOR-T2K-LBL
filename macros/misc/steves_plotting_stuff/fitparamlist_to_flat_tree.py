#!/usr/bin/python

# Turns a TTree containing FitParamLists into an equivalent tree with the FPLs turned into flat variables
# Keeps other tree variables unmodified
# Steve Dennis, 

import ROOT
import math
import array
import numpy as np
import sys
import os

def main(files, treename=None,experiment="NONE"):
  
  try:
    ROOT.valor
  except:
    if   experiment=="T2K" : ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C")
    elif experiment=="DUNE": ROOT.gROOT.ProcessLine(".x $VALOR_DUNE/root/valor_dune.C")
    else:                    ROOT.gROOT.ProcessLine(".x $VALOR_FMWK/root/valor_interpreter.C");
    
  valor = ROOT.valor
  util  = ROOT.valor.util
  conv  = ROOT.valor.conv
  opt   = ROOT.valor.Opt.Instance()
  nr    = ROOT.valor.util.NameRegistry.Instance()

  ext = ".root"

  # Method to call, variable for out_tree, write only if function true on first FPL
  fitparamlist_methods_and_entries= [("BestFitValue","bestfit"  ,""),
                                     ("HesseError"  ,"hesse"    ,"HesseErrorsSet"),
                                     ("MigradError" ,"migraderr","MigradErrorsSet"),
                                     ("InitValue"   ,"init"     ,""),
                                     ("PlusError"   ,"pluserr"  ,"MinosErrorsSet"),
                                     ("MinusError"  ,"minuserr" ,"MinosErrorsSet"),
  ]

  physparamlist_methods_and_entries = [
    ("Value" , "value")
  ]

  # Should be updated to read from Opt::Persistable
  opt.gOscMixingParamChoice = conv.kNuMixSinsqTheta

# REAL MAIN FUNCTION STARTS HERE
  for filename in files:
    
    if not filename.endswith(ext):
      print "Skipping file: ",filename," only supports ",ext," files"
      
    in_file = ROOT.TFile(filename)
    
    base_filename = filename.rsplit("/")[-1]
    output_filename = base_filename[:-len(ext)]+"_flattened.root"
    out_file = ROOT.TFile(output_filename,"RECREATE")
    
    trees = []
    for key in in_file.GetListOfKeys():
      o = key.ReadObj()
      try:
        o.GetListOfBranches()
        trees.append(key.GetName())
      except AttributeError:
        pass
    
    for treename in trees:
      print "Flattening tree ",treename
      in_tree = in_file.Get(treename)
      if not in_tree:
        print "Failed to read tree ",treename," from ",in_file," - skipping"
        continue
      
      out_file.cd()
      
      # Clone the tree without the FitParamLists
      for branch in in_tree.GetListOfBranches():
        if "ParamList" in branch.GetClassName():
          branch.SetStatus(0)
      
      out_tree = in_tree.CloneTree(0)
      
      # Turn the branches back on
      for branch in in_tree.GetListOfBranches():
        if "ParamList" in branch.GetClassName():
          branch.SetStatus(1)
      
      # Get the first tree entry to set up the branches
      in_tree.GetEntry(0)
      
      # Now set up the variables in the flat tree
      new_branches  = []      # Names of the new branches
      var_dict      = dict()  # Branch addresses of new branches
      param_dict    = dict()  # PhysParamId_t of new branches
      basename_dict = dict()  # Name of branch with FitParamList object
      method_dict   = dict()  # Method which must be called on fitparamlist to fill new branch
      
      physics_params = []
      
      # Do FitParamLists
      for branch in in_tree.GetListOfBranches():
        basename = branch.GetName()
        if "FitParamList" in branch.GetClassName():
          fpl = getattr(in_tree,basename) # Get the fitparamlist object
          
          for ipar in xrange(fpl.Size()): # Loop over parameters
            par = fpl.WhichParam(ipar)
            parname = nr.ParamAsStringSimple(par).replace("|","")
            
            if par < conv.kNuis000:
              physics_params.append(par)
            
            for method, entryname,condition in fitparamlist_methods_and_entries:
              if (condition and not (getattr(fpl,condition)())):
                continue
                
              full_branch_name = basename+"__"+parname+"__"+entryname
              new_branches.append(full_branch_name)
              
              var = array.array('d',[0.])
              
              param_dict    [full_branch_name] = par
              var_dict      [full_branch_name] = var
              method_dict   [full_branch_name] = method
              basename_dict [full_branch_name] = basename
              
              out_tree.Branch(full_branch_name,var,full_branch_name+"/D")
      
      # Remove duplicates from physics_params list
      physics_params = list(set(physics_params))
      
      # Do PhysParamLists
      for branch in in_tree.GetListOfBranches():
        basename = branch.GetName()
        if "PhysParamList" in branch.GetClassName():
          ppl = getattr(in_tree,basename) # Get the physparamlist object
          
          # Loop over physics parameters
          for par in xrange(conv.knppmax):
            if not ppl.Contains(par): continue
            parname = nr.ParamAsStringSimple(par).replace("|","")
            
            for method, entryname in physparamlist_methods_and_entries:
              full_branch_name = basename+"__"+parname+"__"+entryname
              new_branches.append(full_branch_name)
              
              var = array.array('d',[0.])
              
              param_dict    [full_branch_name] = par
              var_dict      [full_branch_name] = var
              method_dict   [full_branch_name] = method
              basename_dict [full_branch_name] = basename
              
              out_tree.Branch(full_branch_name,var,full_branch_name+"/D")
              
      # Fill the new tree
      for entry in in_tree:
        for nb in new_branches:
          fpl = getattr(in_tree,basename_dict[nb])                       # Get the required fitparamlist object
          var_dict[nb][0] = getattr(fpl,method_dict[nb])(param_dict[nb]) # Call the method and fill the branch address
          
        out_tree.Fill()
      
      out_tree.Write()
  
if __name__=="__main__":
  # Need argument parser here for experiment - for now check if dune defined
  if "VALOR_DUNE" in os.environ:
    experiment = "DUNE"
  else:
    experiment = "T2K"
    
  main(sys.argv[1:],None,experiment)
