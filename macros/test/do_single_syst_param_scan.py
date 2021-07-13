#!/usr/bin/python

import ROOT

R = ROOT

R.PyConfig.IgnoreCommandLineOptions = True

try:
  R.valor
except:
  R.gROOT.ProcessLine(".x $VALOR_FMWK/root/valor_interpreter.C")
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

#______________________________________________________________________________________________________________________________

def main(config_file, param_values,
    hierarchy,angle_scheme,octant,dmsqscheme,
    hypname,fplopt,
    force_untuned):

  selected_param_name     = "MAQE"
  selected_reac_type      = conv.kmtCCQE
  selected_sample_name    = "superk__JPARC_FHC__1rmu"
  selected_osc_channels   = [
    [conv.kNuMu   , conv.kNuMu],
    [conv.kNuMuBar, conv.kNuMuBar],
    [conv.kNuE    , conv.kNuE],
    [conv.kNuEBar , conv.kNuEBar],
    [conv.kNuMu   , conv.kNuE],
    [conv.kNuMuBar, conv.kNuEBar],
  ]
  low = 0.71
  high = 1.71
  step = 0.01

  opt.gMassHierarchyAssumed = util.HierarchyFromString(hierarchy)
  opt.gOscMixingParamChoice = util.MixingTypeFromString(angle_scheme)
  opt.gOscOctantChoice = util.OctantFromString(octant)
  opt.gOscDmsqScheme = util.DmsqSchemeFromString(dmsqscheme)
  
  opt.SetDefaultParameterValues()

  fit_elements = valor.CompoundFitElements()
  fit_elements.Read(util.CompleteConfigPath(config_file))

  model = fit_elements.FitModel()
  fpl = util.GetFitParameterListFromHypothesisName(hypname, fplopt)

  if force_untuned:
    syst = model.GetSystematics()
    for p in range(0, conv.knppmax):
      if p == conv.kPhypUndefined:
          continue
      if fpl.Exists(p) and syst.Contains(p):
          fpl.Update(p, syst.UntunedValue(p))

  relevant_params = []
  for (parname,val) in param_values.items():
    p = util.NameRegistry.Instance().PhysParamFromString(parname)
    if (fpl.Exists(p)):
      fpl.Update(p,val)
      fpl.SetInitValue(p,val)
      print(util.NameRegistry.Instance().ParamAsString(p).ljust(25),fpl.CurrValue(p))
      relevant_params.append(p)
  relevant_params.sort()



  # Step through linearly spaced values of the selected param and get the event rate for the selected reac type at each step
  selected_param = util.NameRegistry.Instance().PhysParamFromString(selected_param_name)
  vec_val  = ROOT.vector('double')()
  vec_rate = ROOT.vector('double')()
  val = low
  while val <= high:
    fpl.Update(selected_param, val)
    model.Generate(fpl)

    n_submodels = model.GetNSubModels()
    for imodel in range(n_submodels):
      submodel = model.GetSubModel(imodel)
      submodel_name = submodel.Name()
      if submodel_name != selected_sample_name:
        continue

      tag = "valor"
      model_arr = R.TObjArray()
      submodel.Write(model_arr,tag,"pdft_vf pdfc_vf")

      total_for_reac_type = 0
      for osc_channel in selected_osc_channels:
        init_flav  = osc_channel[0]
        final_flav = osc_channel[1]
        m = util.ReacTypeToMode(selected_reac_type, init_flav, final_flav)
        mode_name = util.ModeAsString(m)
        mode_array_name = ("%s_pdfc_ob_vf_%s_%s")%(tag,submodel_name,mode_name)
        spectrum = model_arr.FindObject(mode_array_name)

        assert(spectrum)

        nev = spectrum.Integral()
        total_for_reac_type += nev
      print selected_param_name, "val =", val, ";", util.ReacTypeAsString(selected_reac_type), "rate =", total_for_reac_type
      vec_val.push_back(val)
      vec_rate.push_back(total_for_reac_type)
    val += step

  # Draw
  graph = ROOT.TGraph(vec_val.size(), vec_val.data(), vec_rate.data())
  graph.GetXaxis().SetTitle(selected_param_name)
  graph.GetYaxis().SetTitle(util.ReacTypeAsString(selected_reac_type) + " rate")
  canv = ROOT.TCanvas("canv", "canv", 800, 600)
  graph.Draw("APL")
  canv.SaveAs("single_param_scan_"+selected_param_name+"_"+util.ReacTypeAsString(selected_reac_type)+".pdf")




if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of validation spectra")


  parser.add_option("--physics-hypothesis"      , dest="hypname"       , default="3f_osc"    , help="Name of FPL.")
  parser.add_option("--fit-param-opt"           , dest="fplopt"        , default="NONE"      , help="FPL options.")
  parser.add_option(     '--hierarchy'          , dest="hierarchy"     , default="normal"    , help='Mass hierarchy: "normal" or "inverted"')
  parser.add_option('-m','--mixing-param-choice', dest="mixparamchoice", default="sinsqtheta", help='Mixing param coice: "sinsqtheta" or "sinsq2theta"')
  parser.add_option('-t','--mixing-param-octant', dest="mixparamoctant", default="two"       , help='Mixing param octant: "one" or "two"')
  parser.add_option('-a','--dmsq-scheme'        , dest="dmsqscheme"    , default="auto"      , help='"auto", "abs23", "abs31" or "fogli"')
  
  parser.add_option('-c', '--config', dest='config', default="default_configuration.xml", help="Configuration file")

  parser.add_option('-u', '--force-untuned' , dest='force_untuned' , default=False, action="store_true", help="If true, force all systematics to their untuned values (so their weights should = 1).")

  all_parnames = []
  for p in range(conv.kNuis000):
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
    o.hypname,o.fplopt, o.force_untuned)
  
  if other:
    print("Warning - ignored arguments: ",other)
