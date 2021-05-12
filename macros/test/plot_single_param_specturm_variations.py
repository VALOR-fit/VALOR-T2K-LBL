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

def main(tag, config_file, param_values,
    hierarchy,angle_scheme,octant,dmsqscheme,
    hypname,fplopt,
    force_untuned):

  selected_param_name     = "EB_dial_O_nu"
  single_param_values     = [-10.0    , 0.0       , 10.0     ]
  colours                 = [856      , 1         , 633      ]
  styles                  = [2        , 1         , 2        ]
  labels                  = ["-10 MeV", "No Shift", "+10 MeV"]
  selected_reac_type      = conv.kmtCCQE  # Set to conv.kmtUndefined to look at the sum over all reac types
  selected_kine_var       = conv.kvEreco
  selected_osc_channels   = [
    [conv.kNuMu   , conv.kNuMu],
    [conv.kNuMuBar, conv.kNuMuBar],
    [conv.kNuE    , conv.kNuE],
    [conv.kNuEBar , conv.kNuEBar],
    [conv.kNuMu   , conv.kNuE],
    [conv.kNuMuBar, conv.kNuEBar],
  ]
  
  reac_type_name = util.ReacTypeAsString(selected_reac_type)
  if selected_reac_type == conv.kmtUndefined:
    reac_type_name = "all"
  save_name = "single_param_variation_" +selected_param_name +"_" +reac_type_name +"__" +tag +".pdf"

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

  print "Before stepping through the selected param, the FPL is:"
  fpl.Print()

  # Step through values of the selected param and make the total spectra for the selected reac type at each value
  selected_param = util.NameRegistry.Instance().PhysParamFromString(selected_param_name)
  total_spectra = {}
  for ival, val in enumerate(single_param_values):
    total_spectra[val] = {}
    fpl.Update(selected_param, val)
    model.Generate(fpl,"+SB")

    n_submodels = model.GetNSubModels()
    for imodel in range(n_submodels):
      submodel = model.GetSubModel(imodel)
      submodel_name = submodel.Name()
      total_spectra[val][submodel_name] = None

      sskb = submodel.GetBinning()
      reco_smp = sskb.RecoSmpBinning()
      ks = util.KineSpaceFromKineVars(reco_smp.KineVars())
      ndim = util.NDimReco(ks)

      tag = "valor"
      model_arr = R.TObjArray()
      submodel.Write(model_arr,tag,"pdfc_rf pdfcs_rf pdft_rf pdfts_rf +SB")


      array_name = None
      if selected_reac_type == conv.kmtUndefined:
        if ndim == 1:
          array_name = ("%s_pdft_sb_rf_%s_%s")%(tag,submodel_name,util.KineVarAsString(selected_kine_var))
        else:
          array_name = ("%s_pdfts_sb_rf_%s_%s")%(tag,submodel_name,util.KineVarAsString(selected_kine_var))

        spectrum = model_arr.FindObject(array_name)
        if not spectrum:
          print "ERROR: Could not find spectrum '{}'".format(array_name)
          util.Die(-1)

        if not total_spectra[val][submodel_name]:
          total_spectra[val][submodel_name] = ROOT.TH1D(spectrum)
        else:
          total_spectra[val][submodel_name].Add(spectrum)
        print selected_param_name, "val =", val, ";", submodel_name, "; total spectrum integral =", \
              total_spectra[val][submodel_name].Integral()

      else:
        for osc_channel in selected_osc_channels:
          init_flav  = osc_channel[0]
          final_flav = osc_channel[1]
          m = util.ReacTypeToMode(selected_reac_type, init_flav, final_flav)
          mode_name = util.ModeAsString(m)
          if ndim == 1:
            array_name = ("%s_pdfc_sb_rf_%s_%s_%s")%(tag,submodel_name,util.KineVarAsString(selected_kine_var),mode_name)
          else:
            array_name = ("%s_pdfcs_sb_rf_%s_%s_%s")%(tag,submodel_name,util.KineVarAsString(selected_kine_var),mode_name)

          spectrum = model_arr.FindObject(array_name)
          if not spectrum:
            print "ERROR: Could not find spectrum '{}'".format(array_name)
            util.Die(-1)

          if not total_spectra[val][submodel_name]:
            total_spectra[val][submodel_name] = ROOT.TH1D(spectrum)
          else:
            total_spectra[val][submodel_name].Add(spectrum)
          print selected_param_name, "val =", val, ";", submodel_name, "; channel =", util.FlavourAsString(init_flav), \
                "-->", util.FlavourAsString(final_flav), "; total spectrum integral =", \
                total_spectra[val][submodel_name].Integral()

  # Draw
  n_submodels = model.GetNSubModels()
  for imodel in range(n_submodels):
    submodel = model.GetSubModel(imodel)
    submodel_name = submodel.Name()

    sskb = submodel.GetBinning()
    reco_smp = sskb.RecoSmpBinning()
    ks = util.KineSpaceFromKineVars(reco_smp.KineVars())
    ndim = util.NDimReco(ks)

    canv = ROOT.TCanvas("canv_"+submodel_name, "canv_"+submodel_name, 1000, 1000)

    xmax = 3.0
    if ndim == 2:
      xmax = 1.25

    ymax = 0.0
    for ival, val in enumerate(single_param_values):
      if "1rmu" in submodel_name:
        if "FHC" in submodel_name:
          total_spectra[val][submodel_name].Rebin(2)
        else:
          total_spectra[val][submodel_name].Rebin(3)
      ymax = max(ymax, total_spectra[val][submodel_name].GetMaximum())

    leg_xmin = 0.20
    leg_xmax = 0.35
    if ndim == 1:
      leg_xmin = 0.65
      leg_xmax = 0.80
    leg = ROOT.TLegend(leg_xmin, 0.73, leg_xmax, 0.90)
    leg.SetFillStyle(1001)
    leg.SetBorderSize(0)
    leg.SetTextFont(132)
    leg.SetTextSize(0.04)

    draw_opt = "hist"
    for ival, val in enumerate(single_param_values):
      # Format spectrum
      reac_type_str = reac_type_name
      if selected_reac_type == conv.kmtUndefined:
        reac_type_str = ""
      total_spectra[val][submodel_name].GetXaxis().SetTitle("Reconstructed Neutrino Energy (GeV)")
      total_spectra[val][submodel_name].GetYaxis().SetTitle("Number of " +reac_type_str +" Events")
      #total_spectra[val][submodel_name].GetYaxis().SetTitleOffset(0.9)
      total_spectra[val][submodel_name].GetYaxis().SetTitleSize(0.05)
      total_spectra[val][submodel_name].GetYaxis().SetLabelSize(0.04)
      total_spectra[val][submodel_name].GetXaxis().SetTitleSize(0.05)
      total_spectra[val][submodel_name].GetXaxis().SetLabelSize(0.04)
      total_spectra[val][submodel_name].SetLineColor(colours[ival])
      total_spectra[val][submodel_name].SetLineStyle(styles [ival])
      leg.AddEntry(total_spectra[val][submodel_name], labels[ival], "l")
      total_spectra[val][submodel_name].Draw(draw_opt)
      total_spectra[val][submodel_name].GetXaxis().SetRangeUser(0.0, xmax)
      total_spectra[val][submodel_name].GetYaxis().SetRangeUser(0.0, 1.05*ymax)
      if ival == 0:
        draw_opt += ",same"
    leg.Draw()
    canv.Update()

    if imodel == 0:
      save_name_use = save_name + "("
    elif imodel == n_submodels-1:
      save_name_use = save_name + ")"
    else:
      save_name_use = save_name
    canv.SaveAs(save_name_use)


if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of validation spectra")


  parser.add_option("--tag"                     , dest="tag"           , default="postfit_AsimovA", help="Tag for naming output file.")
  parser.add_option("--physics-hypothesis"      , dest="hypname"       , default="3f_osc"         , help="Name of FPL.")
  parser.add_option("--fit-param-opt"           , dest="fplopt"        , default="NONE"           , help="FPL options.")
  parser.add_option(     '--hierarchy'          , dest="hierarchy"     , default="normal"         , help='Mass hierarchy: "normal" or "inverted"')
  parser.add_option('-m','--mixing-param-choice', dest="mixparamchoice", default="sinsqtheta"     , help='Mixing param coice: "sinsqtheta" or "sinsq2theta"')
  parser.add_option('-t','--mixing-param-octant', dest="mixparamoctant", default="two"            , help='Mixing param octant: "one" or "two"')
  parser.add_option('-a','--dmsq-scheme'        , dest="dmsqscheme"    , default="auto"           , help='"auto", "abs23", "abs31" or "fogli"')
  
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

  main(o.tag, o.config, param_values,
    o.hierarchy,o.mixparamchoice,o.mixparamoctant,o.dmsqscheme,
    o.hypname,o.fplopt, o.force_untuned)
  
  if other:
    print("Warning - ignored arguments: ",other)
