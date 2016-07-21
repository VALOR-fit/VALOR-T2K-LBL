#!/usr/bin/env python

import ROOT

for f in ["sk13a_nRnueCC.root","sk13a_nRnumuCC.root"]:
#for f in ["sk13a_nRnueCC_erec_delta.root","sk13a_nRnumuCC_erec_delta.root"]:
  tf = ROOT.TFile(f)
  
  outfile_2R     = ROOT.TFile(f.replace("nR","2R"    ),"RECREATE")
  outfile_2Rnmue     = ROOT.TFile(f.replace("nR","2Rnmue"    ),"RECREATE")
  outfile_3plusR = ROOT.TFile(f.replace("nR","3plusR"),"RECREATE")
  
  key_names  = [k.GetName() for k in tf.GetListOfKeys()]
  hist_names = [k for k in key_names if k[:2]=="h_"]
  
  for tree_name in ["tnumu","tnue","tnumubar","tnuebar","tnuesig"]:
    tree = tf.Get(tree_name)
    hname = [n for n in hist_names if tree_name in n][0]
    h = ROOT.TH1F(tf.Get(hname))
    
    outfile_2R.cd()
    tree_2R = tree.CloneTree(0)
    outfile_2Rnmue.cd()
    tree_2Rnmue = tree.CloneTree(0)
    outfile_3plusR.cd()
    tree_3plusR = tree.CloneTree(0)
    
    for event in tree:
      if event.nring==2:
        #if event.nmue<=2:
        #tree_2Rnmue.Fill()
        #else:
        tree_2R.Fill()
      else:
        tree_3plusR.Fill()
    
    outfile_2R.cd()
    tree_2R.Write(tree_name)
    h.Write(hname)

    outfile_2Rnmue.cd()
    tree_2Rnmue.Write(tree_name)
    h.Write(hname)
    
    outfile_3plusR.cd()
    outfile_3plusR.Write(tree_name)
    h.Write(hname)
  
  outfile_2R.Close()
  outfile_2Rnmue.Close()
  outfile_3plusR.Close()
  tf.Close()
  
