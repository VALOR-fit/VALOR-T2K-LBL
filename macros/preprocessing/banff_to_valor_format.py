#!/usr/bin/env python
import sys
import ROOT

# Don't let ROOT ruin --help
ROOT.PyConfig.IgnoreCommandLineOptions = True

# load valor
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

# List strange untuned central values. Will assume 1 for all others.
# Want to grab these from the TVector of untuned values instead
abnormal_untuned_central_values = dict()
abnormal_untuned_central_values["CC_DIS"] = 0.0
#abnormal_untuned_central_values["MAQE"] = 1.2
#abnormal_untuned_central_values["CA5"] = 1.01
#abnormal_untuned_central_values["ISO_BKG"] = 1.30
abnormal_untuned_central_values["2p2h_shape_O"] = 0.0
abnormal_untuned_central_values["FEFQE"] = 0.0
abnormal_untuned_central_values["FEFQEH"] = 0.0
abnormal_untuned_central_values["FEFINEL"] = 0.0
abnormal_untuned_central_values["FEFABS"] = 0.0
abnormal_untuned_central_values["FEFCX"] = 0.0
abnormal_untuned_central_values["ebo"] = 0.0

# Set up files
assert(len(sys.argv)==2)

#for beam_mode_name in ["fhc","rhc","all"]:
for beam_mode_name in ["all"]:
  infile_name = sys.argv[1]
  
  outfile_name = infile_name.replace(".root","_valor_order_"+beam_mode_name+".root")
  print "File names", infile_name, outfile_name

  infile  = ROOT.TFile(infile_name)
  outfile = ROOT.TFile(outfile_name,"RECREATE")

  # Copy the things we don't have to change
  for key in infile.GetListOfKeys():
    key_name = key.GetName()
    if key_name not in ["param_list","prefit_params","prefit_cov","postfit_params","postfit_cov"]:
      print "Copying", key_name
      o = key.ReadObj()
      outfile.cd()
      o.Write(key_name)
      
  # Get the input parameters and count out the ones we need
  input_param_list = [str(s) for s in infile.Get("param_list")]
  n_banff_flux_params = 0
  for s in input_param_list:
    if "SKNuMode" in s or "SKANuMode" in s:
      n_banff_flux_params += 1

  n_banff_params = len(input_param_list)
  n_banff_xsec_params = sum([int(conv.kXSecSystIncludedInBANFF[i]) for i in xrange(conv.kNumOfXSecSyst)])
  
  # REMOVE SCC PARAMETERS FROM XSEC - FOR inject_scc VERSION ONLY
  for position_in_valor_list in xrange(0,conv.kNumOfFHCpRHCParamsInBANFFfit): # just use the FHC versions here, will need to change if RHC/FHC systematics differ
    p = util.GetBANFFParam(position_in_valor_list)
    if "banffscc" in util.ParamAsStringCompact(p): n_banff_xsec_params -= 1

  # Hack to get around additional parameters in the BANFF matrix
  # which aren't in the analysis

  first_xsec_param = n_banff_params - ( n_banff_xsec_params )
  
  # n_banff_flux_params = n_banff_params - ( n_banff_xsec_params + additional_xsec_params_in_matrix )

  print "n_banff_params", n_banff_params  
  print "n_banff_flux_params", n_banff_flux_params
  print "n_banff_xsec_params", n_banff_xsec_params, "from", first_xsec_param

  
  assert(n_banff_flux_params%2 == 0)
  
  n_banff_flux_params_per_mode = n_banff_flux_params/2
  if (beam_mode_name == "fhc"):
    first_flux_param = 0
    last_flux_param  = n_banff_flux_params_per_mode-1
  if (beam_mode_name == "rhc"):
    first_flux_param = n_banff_flux_params_per_mode
    last_flux_param  = first_xsec_param-1
  if (beam_mode_name == "all"):
    first_flux_param = 0
    last_flux_param  = first_xsec_param-1

  print "All parameters"
  print input_param_list
  print "Flux parameters"
  print input_param_list[first_flux_param:last_flux_param+1]
  print "xsec parameters"
  print input_param_list[first_xsec_param:n_banff_params]
  
  # Build a list of the order of the parameters we want to select
  # first grab all the flux params, in order
  param_order_list = list(xrange(first_flux_param,last_flux_param+1))
  print "Param ordered list", param_order_list
  
  # Now look at xsec params
  # just use the FHC versions here, will need to change if RHC/FHC systematics differ
  for position_in_valor_list in xrange(0,conv.kNumOfFHCpRHCParamsInBANFFfit): 
    phys_param_id = util.GetBANFFParam(position_in_valor_list)
    param_name = util.ParamAsStringNIWG(phys_param_id)
    if param_name == "--": continue
    if "SCC" in param_name: continue
    position_in_banff_list = input_param_list.index(param_name)
    param_order_list.append(position_in_banff_list)
    print "iv", phys_param_id, "V:", util.ParamAsString(phys_param_id), "N:", param_name, "inpos:", position_in_banff_list

  print "Complete ordered list", param_order_list

  reordered_list = [ input_param_list[i] for i in param_order_list]

  print "Param Order List after XSec params added (%s)" % beam_mode_name
  print reordered_list
  output_names = ROOT.TObjArray()
  for p in param_order_list:
    output_names.Add(ROOT.TObjString(input_param_list[p]))
    
  output_names.Write("param_list",ROOT.TObject.kSingleKey)
  
  n_elements = len(param_order_list) 
  print "Elements in new matrix", n_elements
  
  reordered_covs_d = dict()
  reordered_params_d = dict()
  
  # Reorder the parameters and the covariance
  for fit in ["prefit","postfit"]:
    orig_cov    = infile.Get(fit+"_cov")
    orig_params = infile.Get(fit+"_params")

    reordered_cov    = ROOT.TMatrixD(n_elements,n_elements)
    reordered_params = ROOT.TVectorD(n_elements)
    for i,orig_i in enumerate(param_order_list):
      reordered_params[i] = orig_params[orig_i]
      for j,orig_j in enumerate(param_order_list):
        reordered_cov[i][j] = orig_cov[orig_i][orig_j]
        
    outfile.cd()
    reordered_cov   .Write(fit+"_cov")
    reordered_params.Write(fit+"_params")
    
    reordered_covs_d  [fit] = reordered_cov
    reordered_params_d[fit] = reordered_params
    
    
  # Produce untuned central values - all
  vec = ROOT.TVectorD(n_banff_params)
  
  for i in xrange(vec.GetNrows()):
    vec[i] = 1.
  for name, val in abnormal_untuned_central_values.iteritems():
    try:
      ind = list(output_names).index(name)
      vec[ind] = val
    except ValueError:
      print "Making central values - failed to set value of ",name
      
  vec.Write("banff_untuned_central_values")
  
  
  # Produce untuned central values - xsec only
  xsec_vec = ROOT.TVectorD(n_banff_xsec_params)
  xsec_output_names = ROOT.TObjArray(output_names)
  while xsec_output_names.GetEntries() > n_banff_xsec_params:
    xsec_output_names.RemoveFirst()
    xsec_output_names.Compress()
  xsec_output_names.Write("param_list_banff_xsec",ROOT.TObject.kSingleKey)
  
  for i in xrange(xsec_vec.GetNrows()):
    #xsec_vec[i] = 
    xsec_vec[i] = 1.
  for name, val in abnormal_untuned_central_values.iteritems():
    try:
      ind = list(xsec_output_names).index(name)
      xsec_vec[ind] = val
    except ValueError:
      print "Making central values - failed to set value of ",name
  
  xsec_vec.Write("banff_xsec_untuned_central_values")
    
  infile.Close()
  outfile.Close()
  
  # Start SCC injection
  injected_outfile_name = outfile_name.replace(".root","_plus_scc.root")
  print "Writing output file", injected_outfile_name
  injected_outfile = ROOT.TFile(injected_outfile_name,"RECREATE")
  injected_output_names = ROOT.TObjArray()
  for i,p in enumerate(param_order_list):
    injected_output_names.Add(ROOT.TObjString(input_param_list[p]))
    if input_param_list[p] == "MARES":
      injected_output_names.Add(ROOT.TObjString("SCCA"))
      injected_output_names.Add(ROOT.TObjString("SCCV"))
      injection_point = i
  injected_output_names.Write("param_list",ROOT.TObject.kSingleKey)
  
  scc_val = 1
  
  n_elements_i = n_elements+2
  for fit in ["prefit","postfit"]:
    reordered_params = reordered_params_d[fit]
    reordered_cov    = reordered_covs_d  [fit]
    
    injected_cov    = ROOT.TMatrixD(n_elements_i,n_elements_i)
    injected_params = ROOT.TVectorD(n_elements_i)
    i = 0
    i_old = 0
    while i < n_elements_i:
      injected_params[i] = reordered_params[i_old]
      j = 0
      j_old = 0
      while j < n_elements_i:
        injected_cov[i][j] = reordered_cov[i_old][j_old]
        if j==injection_point:
          j+=2
        j    +=1
        j_old+=1
        
      if i==injection_point:
        injected_params[i+1] = scc_val
        injected_params[i+2] = scc_val
        i+=2
      i+=1
      i_old+=1
      
    injected_cov[injection_point+1][injection_point+1] = scc_val
    injected_cov[injection_point+2][injection_point+2] = scc_val
    
    injected_outfile.cd()
    injected_cov   .Write(fit+"_cov")
    injected_params.Write(fit+"_params")
  
  # Produce injected untuned central values - all
  injected_vec = ROOT.TVectorD(n_banff_params + 2)
  
  for i in xrange(injected_vec.GetNrows()):
    injected_vec[i] = 1.
  for name, val in abnormal_untuned_central_values.iteritems():
    try:
      ind = list(injected_output_names).index(name)
      injected_vec[ind] = val
    except ValueError:
      print "Making SCC injected central values - failed to set value of ",name
  
  injected_vec.Write("banff_untuned_central_values")
  
  # Produce injected untuned central values - xsec only
  injected_xsec_vec = ROOT.TVectorD(n_banff_xsec_params + 2)
  injected_xsec_output_names = ROOT.TObjArray(injected_output_names)
  while injected_xsec_output_names.GetEntries() > n_banff_xsec_params+2:
    injected_xsec_output_names.RemoveFirst()
    injected_xsec_output_names.Compress()
  injected_xsec_output_names.Write("param_list_banff_xsec",ROOT.TObject.kSingleKey)
  
  for i in xrange(injected_xsec_vec.GetNrows()):
    injected_xsec_vec[i] = 1.
  for name, val in abnormal_untuned_central_values.iteritems():
    try:
      ind = list(injected_xsec_output_names).index(name)
      injected_xsec_vec[ind] = val
    except ValueError:
      print "Making SCC injected central values - failed to set value of ",name
  
  injected_xsec_vec.Write("banff_xsec_untuned_central_values")
  
  injected_outfile.Close()
