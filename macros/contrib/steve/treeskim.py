#!/usr/bin/env python

import ROOT as R

R.gROOT.ProcessLine(".x $GENIE/src/scripts/gcint/loadlibs.C");
#~ R.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C");
#~ 
#~ valor = R.valor

def accept_event(ev):
  if ev.sample>0:
    return True
  return False

def main(infiles,outfile_pattern,events_per_file,treename):
  chain = R.TChain(treename)
  for f in infiles:
    chain.Add(f)
  
  assert(chain.GetEntries())
  print "Tree has ",chain.GetEntries," entries."
  
  output_file = None
  output_file_number = 0
  
  out_tree = None
  
  events_written = 0
  
  for i,event in enumerate(chain):
    if (events_written%events_per_file==0):
      outfile_name = outfile_pattern.replace("[T]",treename).replace("[N]",str(output_file_number))
      if output_file is None or "[N]" in outfile_pattern:
        if out_tree and output_file: out_tree.Write()
        if output_file: output_file.Close()
        output_file = R.TFile(outfile_name,"RECREATE")
        output_file_number+=1
        out_tree = chain.CloneTree(0)
    if accept_event(event):
      out_tree.Fill()
      events_written += 1
    if i%1000==0:
      print "Studied",i+1," events. Selected",events_written
  
  out_tree.Write()
    
  print "Done!"

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Submit many VALOR sensitivity jobs")
  
  parser.add_option('-o', '--outfile', dest='outfile', default="skimmed_[T]_[N].root", help="Output ROOT file")
  parser.add_option('-n', '--nevents-per-file', dest='events_per_file', default=100000,type=int, help="Number of events per output ROOT file")
  parser.add_option('-t', '--treename', dest='treename', default="fitnt",help="Name of tree to skim.")
  
  o,args = parser.parse_args()
  
  main(args,o.outfile,o.events_per_file,o.treename)
