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
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()
arridx= util.ArrIdx.Instance()

species        = ["numu","nue","numubar","nuebar","signue", "signuebar"]
processes      = ["CCQE", "CC1pi", "CCcoherent", "CCother", "NC1pi0" , "NC1piPM", "NCcoherent", "NCother"]
processes_used = ["CCQE", "CC1pi", "CCcoherent", "CCother", "NC1pi0" , "NC1piPM", "NCcoherent", "NCother"]
#~ processes_used = ["CCQE", "CC1pi", "CCcoherent", "CCother" , "NC1piPM", "NCother"]

col1_width = 30
coln_width = 3

include_empty = True

#______________________________________________________________________________________________________________________________

def main(sample=conv.kSKsFCFV1RingMuLike,output_file="syst_applied_table.tex", skipflux=False, skipskdet=False):
  
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/analysis_configuration.xml")
  skdetfsi = skexp.GetSuperkDetectorFSI()
  xsec = skexp.GetNuXSec()

  print "Sample enum: ", sample 
  
  # Some ridiculous voodoo to go around ROOT's inability to identify arrays of enums or 2D arrays of bools
  modes_array  = array.array("i",[-2]*conv.knm)
  
  #~ ROOT.gROOT.ProcessLine("memcpy("+str(modes_array.buffer_info()[0])+", &(valor::conv::kSKReactionModes["+`sample`+"][0]),valor::conv::knm*sizeof(int))")
  root_loop = "{using namespace valor::conv; int * arr = "+str(modes_array.buffer_info()[0])
  root_loop += ";for (int i = 0 ; i < valor::conv::knm;i++){arr[i]=valor::conv::kSKReactionModes["+`sample`+"][i];};};"
  ROOT.gROOT.ProcessLine(root_loop)
  
  modes = [m for m in modes_array if m != conv.kSKmUndefined]
  
  #~ soff = sample*conv.knm
  
  tex_output = open(output_file,"w")
  
  tex_output.write("""\documentclass[10pt]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{fullpage}
\\usepackage{rotate}  
\\usepackage{rotating}
\\usepackage[paperwidth=300mm,paperheight=600mm]{geometry}

\\usepackage[table]{xcolor} % colored cells for the output table

\\newcommand{\ck}{\\cellcolor{green!50}$\\checkmark$}
\\begin{document}

{
\\setlength{\\tabcolsep}{3pt}
\\renewcommand{\\arraystretch}{1.8}
""")


  # Write osc params
  phys_params = [conv.kPhypMixing12, conv.kPhypMixing13, conv.kPhypDm12sq, conv.kPhypCPV]
  
  tex_output.write("\\begin{table} \\tiny\n")
  #~ tex_output.write("\\centering\n")
  tex_output.write("\n\n\\begin{tabular} { | c ||"+("c | "*len(modes))+" | }\n\\hline\n")
  
  for mode in modes:
    tex_output.write(" & \\begin{sideways} "+util.ModeAsStringTex(mode)+" \end{sideways} \n")
  
  tex_output.write("\\\\ \hline \n")
  
  # Physics
  for param in phys_params:
    entries = [util.ParamAsStringTex(param).ljust(col1_width)]
    for m in modes:
      if util.IsCC(m):
        entries.append("\\ck".ljust(coln_width))
      else:
        entries.append("".ljust(coln_width))
    tex_output.write("& ".join(entries) + " \\\\ \hline \n")
      
  # SKDetFSI
  if not skipskdet:
    for param in xrange(conv.kNuis000,conv.kLastNuisParam):
      if not util.IsSKDetFSI(param): continue
      entries = [util.ParamAsStringPaperTex(param).ljust(col1_width)]
      row = False
      for m in modes:
        if skdetfsi.ApplySystToSKMode(param,sample,m):
          entries.append("\\ck".ljust(coln_width))
          row = True
        else:
          entries.append("".ljust(coln_width))
      if row or include_empty: tex_output.write("& ".join(entries) + " \\\\ \hline \n")
    
  # BANFF
  for param in xrange(conv.kNuis000,conv.kLastNuisParam):
    if not util.IsBANFF(param): continue
    entries = [util.ParamAsStringPaperTex(param).ljust(col1_width)]
    row = False
    if util.IsRHCSample(sample): ip = arridx.RHCBANFFSystIdx(param)
    else:                        ip = arridx.BANFFSystIdx(param)
    if ip < 0: continue
    
    if conv.kParamInBANFFfitIsSuperkFluxParam[ip]: # flux parameters
      if skipflux: continue
      nupdg = conv.kNuPdgCodeOfFluxParamInBANFFfit[ip]
      assert(nupdg)
      for m in modes:
        if util.SuperkMode2FluxNuPdg(m)==nupdg:
          entries.append("\\ck".ljust(coln_width))
          row = True
        else:
          entries.append("".ljust(coln_width))
      if row or include_empty: tex_output.write("& ".join(entries) + " \\\\ \hline \n")
    
    else: #banff xsec parameters
      row = False
      for m in modes:
        if xsec.ApplySystToSKMode(param,sample,m):
          entries.append("\\ck".ljust(coln_width))
          row = True
        else:
          entries.append("".ljust(coln_width))
      if row or include_empty: tex_output.write("& ".join(entries) + " \\\\ \hline \n")
      
  # XSec
  for param in xrange(conv.kNuis000,conv.kLastNuisParam):
    if not util.IsXSec(param): continue
    
    if util.IsRPA(param):
      print "Skipping RPA param: ",util.ParamAsString(param)
      continue
      
    entries = [util.ParamAsStringPaperTex(param).ljust(col1_width)]
    row = False
    for m in modes:
      if xsec.ApplySystToSKMode(param,sample,m):
        entries.append("\\ck".ljust(coln_width))
        row = True
      else:
        entries.append("".ljust(coln_width))
    if row: tex_output.write("& ".join(entries) + " \\\\ \hline \n")
    if not row: continue
    
  tex_output.write("\n\\end{tabular}\n")
  tex_output.write("\\end{table}\n}\n")
  tex_output.write("\\end{document}\n")
  
  print "\n\nTable output written to ",output_file

if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of applied systematics")
  
  parser.add_option('-f', '--file', dest='output_file', default="syst_applied_table.tex", help="Name for output TeX file.")
  
  parser.add_option('-s', '--sample' , dest='sample'      , default="RHC1Re"                    , help="Sample (1Rmu, 1Re, RHC1Rmu, RHC1Re, MultiRmu, MultiRe)")
  parser.add_option('--skip-skdetfsi', dest='skipskdetfsi', action="store_true", default=False, help="Skip the SKdetFSI parameters")
  parser.add_option('--skip-flux'    , dest='skipflux'    , action="store_true", default=False, help="Skip the flux parameters")
  
  o,other = parser.parse_args()
  sample = util.SKSampleFromString(o.sample)
  main(sample,o.output_file, o.skipflux, o.skipskdetfsi)
  
  if other:
    print "Warning - ignored arguments: ",other
