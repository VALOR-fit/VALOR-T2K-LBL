#!/usr/bin/python

import ROOT
import math
import array
import sys
import os

R = ROOT

R.PyConfig.IgnoreCommandLineOptions = True

try:
  R.valor
except:
  R.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C");
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

col1_width = 12
coln_width = 20
min_channels_to_landscape = 7   # If we get more than this number of channels, then switch the page to landscape orientation

#______________________________________________________________________________________________________________________________

def main(config_file, param_values,
    hierarchy,angle_scheme,octant,dmsqscheme,
    output_file,label,hypname,fplopt,num_decimal_places,
    force_untuned):
  
  opt.gMassHierarchyAssumed = util.HierarchyFromString(hierarchy)
  opt.gOscMixingParamChoice = util.MixingTypeFromString(angle_scheme)
  opt.gOscOctantChoice = util.OctantFromString(octant)
  opt.gOscDmsqScheme = util.DmsqSchemeFromString(dmsqscheme)
  
  opt.SetDefaultParameterValues()

  to_ndp = "{0:."+str(num_decimal_places)+"f}"

  fit_elements = valor.CompoundFitElements()
  fit_elements.Read(util.CompleteConfigPath(config_file))

  model = fit_elements.FitModel()
  fpl = util.GetFitParameterListFromHypothesisName(hypname, fplopt)

  if force_untuned:
    model.ForceSystematicsToUntuned()

  relevant_params = []
  for (parname,val) in param_values.iteritems():
    p = util.NameRegistry.Instance().PhysParamFromString(parname)
    if (fpl.Exists(p)):
      fpl.Update(p,val)
      fpl.SetInitValue(p,val)
      print util.NameRegistry.Instance().ParamAsString(p).ljust(25),fpl.CurrValue(p)
      relevant_params.append(p)
  relevant_params.sort()
  
  tex_output = open(output_file,"w")
  
  tex_output.write("""\documentclass[10pt,a4paper]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{fullpage}
\\usepackage{lscape}
\\usepackage{geometry}
\\geometry{margin=1in}
\\begin{document}
""")

  # Write osc params
  
  tex_output.write("\n\\section*{Job details}\n\n")
  
  tex_output.write("\n\\begin{minipage}[t]{0.5\\textwidth}")
  if len(relevant_params):
    tex_output.write("\n\\begin{itemize}\n")
    for param in relevant_params:
      tex_output.write("\\item[] "+util.NameRegistry.Instance().ParamAsStringTex(param)+" = "+str(fpl.CurrValue(param))+"\n")

    tex_output.write("\\end{itemize}\n")

  tex_output.write("\\end{minipage}")

  tex_output.write("\n\\\\\n\n")
  phyp = fpl.PhysicsHypothesis
  if (phyp=="3f_osc") :
    tex_output.write("Hierarchy: "+hierarchy.capitalize()+"\n\n")
  if label: tex_output.write(label+"\n\n")

  model.Generate(fpl)
  
  n_submodels = model.GetNSubModels()
  
  for imodel in xrange(n_submodels):
    
    submodel = model.GetSubModel(imodel)
    submodel_name = submodel.Name()
    tag = "valor"
    model_arr = R.TObjArray()
    submodel.Write(model_arr,tag,"pdft_vf pdfc_vf")
    
    nev_tot = 0.
    
    modes = [m for m in submodel.GetReactionModes() if m >=0 and m<1E4]
    
    reac_types = []
    channels = [] # here channel is init neutrino flavour + final neutrino flavour
    
    for m in modes:
      init_flav  = util.ReacModeToInitFlavor(m)
      final_flav = util.ReacModeToFinalFlavor(m)
      reac_type  = util.ReacModeToType(m)
      if reac_type not in reac_types: reac_types.append(reac_type)
      channel = (init_flav,final_flav)
      if channel not in channels: channels.append(channel)

    if len(channels) > min_channels_to_landscape:
      tex_output.write("\\begin{landscape}\n")

      
    event_counts_per_channel         = dict() # Holds dictionaries of formatted strings
    total_event_counts_per_channel   = dict() # Holds float
    total_event_counts_per_reac_type = dict() # Holds float
    
    for reac_type in reac_types:
      event_counts_per_channel        [reac_type] = dict()
      total_event_counts_per_reac_type[reac_type] = 0.
      
    for channel in channels:
      total_event_counts_per_channel[channel] = 0.
      
    # Signal not implemented yet
    signal_events     = 0.
    background_events = 0.
    
    for m in modes:
      
      mode_name = util.ModeAsString(m)
      mode_array_name = ("%s_pdfc_ob_vf_%s_%s")%(tag,submodel_name,mode_name)
      spectrum = model_arr.FindObject(mode_array_name)
      
      assert(spectrum)
      
      channel = (util.ReacModeToInitFlavor(m),util.ReacModeToFinalFlavor(m))
      reac_type = util.ReacModeToType(m)
      
      nev = spectrum.Integral()
      nev_tot += nev
      total_event_counts_per_channel  [channel] += nev
      total_event_counts_per_reac_type[reac_type] += nev
      event_counts_per_channel[reac_type][channel] = to_ndp.format(nev)
      
      # ~ if (signal_for_sample[sample_name]==spec): signal_events+=nev
      # ~ else: background_events+=nev
      
      # ~ print channel,"\t",util.ReacTypeAsString(reac_type),reac_type,"\t",event_counts_per_channel[reac_type][channel]
      
    
    # Format the chanenls in a nice way
    channel_pattern = r"$%s \rightarrow %s$"
    fancy_channels = dict([(c,channel_pattern%tuple(map(util.FlavourAsStringTeX,c))) for c in channels])
    fancy_channels["total"] = "Total"
    
    # Start writing output
    submodel_name_escaped_underscores = submodel_name.replace("_",r"\_")
    tex_output.write("\n\\section*{"+submodel_name_escaped_underscores+"}\n")
    tex_output.write("\\begin{tabular} { | l ||"+("r | "*len(channels))+"| r | }\n\\hline\n")
    tex_output.write(" "+"".ljust(col1_width)+" & "+ " & ".join([fancy_channels[c].ljust(coln_width) for c in channels+["total"]]) + " \\\\ \\hline\\hline\n")
    
    for reac_type in reac_types:
      ev = event_counts_per_channel[reac_type]
      # ~ if not total_event_counts_per_reac_type[reac_type]: continue
      reac_type_name = util.ReacTypeAsString(reac_type)
      tex_output.write(" "+reac_type_name.ljust(col1_width)+" & ")
      ev_counts = []
      for channel in channels:
        try:
          ev_counts.append(ev[channel].ljust(coln_width))
        except KeyError:
          ev_counts.append("N/A".ljust(coln_width))
      
      ev_counts.append(to_ndp.format(total_event_counts_per_reac_type[reac_type]).ljust(coln_width))
          
      tex_output.write(" & ".join(ev_counts) + " \\\\ \\hline\n")
    
    # Write totals
    tex_output.write("\\hline\n")
    tex_output.write(" "+"Total".ljust(col1_width)+" & ")
    string_totals = [to_ndp.format(total_event_counts_per_channel[channel]).ljust(coln_width) for channel in channels]
    string_totals.append(to_ndp.format(nev_tot).ljust(coln_width))
    tex_output.write(" & ".join(string_totals) + " \\\\ \\hline\n")
    
    tex_output.write("\n\\end{tabular}\n")
    
    total_array_name = ("%s_pdft_ob_vf_%s")%(tag,submodel_name)
    
    total_spectrum  = model_arr.FindObject(total_array_name)
    assert(total_spectrum)
    total_events = total_spectrum.Integral()
    
    print submodel_name,total_events
    
    tex_output.write("\nTotal: "           +to_ndp.format(total_events)     +"    ("+to_ndp.format(nev_tot)+")\n\n")
    # ~ tex_output.write("\nTotal signal: "    +to_ndp.format(signal_events)    +"\n\n")
    # ~ tex_output.write("\nTotal background: "+to_ndp.format(background_events)+"\n\n")
    
    tex_output.write("POT: "+to_ndp.format(submodel.GetPOT()*1E-20)+"E20\n\n")

    if len(channels) > min_channels_to_landscape:
      tex_output.write("\\end{landscape}\n")

    tex_output.write("\n\n\\clearpage\n\n")
    
  
  tex_output.write("\\end{document}\n")
  
  print "Table output written to ",output_file

if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of validation spectra")


  parser.add_option("--physics-hypothesis"      , dest="hypname"       , default="3f_osc"    , help="Name of FPL.")
  parser.add_option("--fit-param-opt"           , dest="fplopt"        , default="NONE"      , help="FPL options.")
  parser.add_option(     '--hierarchy'          , dest="hierarchy"     , default="normal"    , help='Mass hierarchy: "normal" or "inverted"')
  parser.add_option('-m','--mixing-param-choice', dest="mixparamchoice", default="sinsqtheta", help='Mixing param coice: "sinsqtheta" or "sinsq2theta"')
  parser.add_option('-t','--mixing-param-octant', dest="mixparamoctant", default="two"       , help='Mixing param octant: "one" or "two"')
  parser.add_option('-a','--dmsq-scheme'        , dest="dmsqscheme"    , default="auto"      , help='"auto", "abs23", "abs31" or "fogli"')
  
  parser.add_option('-f', '--file', dest='output_file', default="validation.tex", help="Name for output TeX file.")
  
  parser.add_option('-c', '--config', dest='config', default="default_configuration.xml", help="Configuration file")
  
  parser.add_option('-l', '--label' , dest='label' , default=None, help="Label to put in file.")

  parser.add_option('-d', '--n-decimal-places' , dest='ndp' , default=3, type=int, help="Number of decimal places used to print event rates.")

  parser.add_option('-u', '--force-untuned' , dest='force_untuned' , default=False, action="store_true", help="If true, force all systematics to their untuned values (so their weights should = 1).")

  all_parnames = []
  for p in xrange(conv.kNuis000):
    parname = util.NameRegistry.Instance().ParamAsStringVerySimple(p)
    default = opt.GetGlobalBestFit(p)
    parser.add_option('--%s'%parname, dest=parname,default=default,type=float, help='Value of %s'%parname)
    all_parnames.append(parname)

  v = optparse.Values()
  parser.parse_args(values=v) # v contains only the values that were specified, not defaults
  
  o,other = parser.parse_args()
  
  # Make a dictionary of the parameter names to their argument values
  param_values = dict(((pname,getattr(o,pname)) for pname in all_parnames))

  main(o.config, param_values,
    o.hierarchy,o.mixparamchoice,o.mixparamoctant,o.dmsqscheme,
    o.output_file,o.label,o.hypname,o.fplopt,o.ndp, o.force_untuned)
  
  if other:
    print "Warning - ignored arguments: ",other
