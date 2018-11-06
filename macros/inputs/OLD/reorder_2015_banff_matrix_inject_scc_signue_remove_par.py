#!/usr/bin/env python
import sys
import ROOT


#
# THIS MATRIX WORKS WELL ONLY FOR THE FULL FHC+RHC COV. MATRIX
# FOR SUB-MATRICES THE REMOVAL OF PARAMETERS DOESN'T WORK
# BECAUSE THE INDICES (param_reject_ind) ARE DIFFERENT FOR EACH MATRIX
#

# Set inputs

param_reject_ind = dict()
param_reject_name = dict()
param_reject_ind = [ # original matrix order
  #57 # remove CCNUE_0
]
param_reject_name = [ # original matrix name
  #"CCNUE_0" # remove CCNUE_0
] 



index_new_orig = [] # cov matrix index of valor and original ones

for i in range(50):
  index_new_orig.append([i,i]) # the flux order is identical

# index_new_orig.append([50,52]) # MEC_O
# index_new_orig.append([51,54]) # CA5
# index_new_orig.append([52,56]) # BgRES
# index_new_orig.append([53,50]) # MAQE
# index_new_orig.append([54,55]) # MANFFRES
# index_new_orig.append([55,-99]) # SCCA
# index_new_orig.append([56,-99]) # SCCV
# index_new_orig.append([57,51]) # pF_O
# index_new_orig.append([58,58]) # DISMPISHP
# index_new_orig.append([59,53]) # EB_O
# index_new_orig.append([60,59]) # CCCOH_O_0
# index_new_orig.append([61,60]) # NCCOH_0
# index_new_orig.append([62,62]) # NCOTHER_FAR_0
# index_new_orig.append([63,64]) # CCNue
# index_new_orig.append([64,61]) # NC1GAMMA_0
# index_new_orig.append([65,65]) # CCNueBar
# index_new_orig.append([66,63]) # MEC_NUBAR

index_new_orig.append([50,52]) # MEC_O
index_new_orig.append([51,54]) # CA5
index_new_orig.append([52,56]) # BgRES
index_new_orig.append([53,50]) # MAQE
index_new_orig.append([54,55]) # MANFFRES
index_new_orig.append([55,-99]) # SCCA
index_new_orig.append([56,-99]) # SCCV
index_new_orig.append([57,51]) # pF_O
index_new_orig.append([58,57]) # DISMPISHP
index_new_orig.append([59,53]) # EB_O
index_new_orig.append([60,58]) # CCCOH_O_0
index_new_orig.append([61,59]) # NCCOH_0
index_new_orig.append([62,61]) # NCOTHER_FAR_0
index_new_orig.append([63,63]) # CCNue
index_new_orig.append([64,60]) # NC1GAMMA_0
index_new_orig.append([65,64]) # CCNueBar
index_new_orig.append([66,62]) # MEC_NUBAR



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
abnormal_untuned_central_values = dict()
abnormal_untuned_central_values["DISMPISHP"] = 0.

# Set up files
assert(len(sys.argv)==2)

#for beam_mode_name in ["fhc","rhc","all"]:
for beam_mode_name in ["all"]:
  infile_name = sys.argv[1]
  
  outfile_name = infile_name.replace(".root","_valor_order_"+beam_mode_name+".root")

  infile  = ROOT.TFile(infile_name)
  outfile = ROOT.TFile(outfile_name,"RECREATE")

  # Copy the things we don't have to change
  for key in infile.GetListOfKeys():
    key_name = key.GetName()
    if key_name not in ["param_list","prefit_params","prefit_cov","postfit_params","postfit_cov"]:
      o = key.ReadObj()
      outfile.cd()
      o.Write(key_name)
      




  # # Get the input parameters and count out the ones we need
  # input_param_list = [str(s) for s in infile.Get("param_list")]

  # NEW
  dict_input_param_list = dict();
  input_param_list = []
  i=0
  for s in infile.Get("param_list"):

    #if s=="CCNUE_0": continue    
    if s in param_reject_name:
      continue

    dict_input_param_list[i] = str(s)
    input_param_list.append(str(s))
    i+=1

  #for j in xrange(0,len(input_param_list)):
  #  print j," --> ",input_param_list[j]
  #print input_param_list.index("MAQE")
  #exit(1)
  





  n_banff_params = len(input_param_list)
  n_banff_xsec_params = sum([int(conv.kXSecSystIncludedInBANFF[i]) for i in xrange(conv.kNumOfXSecSyst)])

  # REMOVE SCC PARAMETERS FROM XSEC - FOR inject_scc VERSION ONLY
  for position_in_valor_list in xrange(0,conv.kNumOfFHCpRHCParamsInBANFFfit): # just use the FHC versions here, will need to change if RHC/FHC systematics differ
    p = util.GetBANFFParam(position_in_valor_list)
    if "banffscc" in util.ParamAsStringCompact(p): n_banff_xsec_params -= 1
   
  first_xsec_param = n_banff_params-n_banff_xsec_params
  
  n_banff_flux_params = n_banff_params - n_banff_xsec_params
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
  
  # Build a list of the order of the parameters we want to select
  # first grab all the flux params, in order
  param_order_list = list(xrange(first_flux_param,last_flux_param+1))
  
  # Now look at xsec params
  for position_in_valor_list in xrange(0,conv.kNumOfFHCpRHCParamsInBANFFfit): # just use the FHC versions here, will need to change if RHC/FHC systematics differ
    phys_param_id = util.GetBANFFParam(position_in_valor_list)
    param_name = util.ParamAsStringNIWG(phys_param_id)
    if param_name == "--": continue
    if "SCC" in param_name: continue

    position_in_banff_list = input_param_list.index(param_name)
    param_order_list.append(position_in_banff_list)
    

  output_names = ROOT.TObjArray()
  for p in param_order_list:
    output_names.Add(ROOT.TObjString(input_param_list[p]))
    
  output_names.Write("param_list",ROOT.TObject.kSingleKey)

  n_elements = len(param_order_list)
  
  reordered_covs_d = dict()
  reordered_params_d = dict()
  
  # Reorder the parameters and the covariance
  for fit in ["prefit","postfit"]:
    orig_cov    = infile.Get(fit+"_cov")
    orig_params = infile.Get(fit+"_params")
    
    reordered_cov    = ROOT.TMatrixD(n_elements,n_elements)
    reordered_params = ROOT.TVectorD(n_elements)
    



    for i,orig_i in enumerate(param_order_list): # param_order_list has VALOR order and # of params
      
      # NEW
      n_wrong_i = 0 
      flag_wrong_i = 0
      for iwrong in param_reject_ind:
        if orig_i >= iwrong:  
          n_wrong_i += 1
      #print "iwrong =",iwrong," - n_wrong_i =",n_wrong_i,"i =",i," - orig_i =",orig_i
      orig_i=orig_i+n_wrong_i # SKIP WRONG PARAMETERS
      #if orig_i>=57:  # SKIP WRONG PARAMETER
      #  orig_i=orig_i+1       
      #

      reordered_params[i] = orig_params[orig_i]

            
      for j,orig_j in enumerate(param_order_list):

        # NEW
        n_wrong_j = 0 
        flag_wrong_j = 0
        for jwrong in param_reject_ind:
          if orig_j >= jwrong:  
            n_wrong_j += 1
        #print "n_wrong_j = ",n_wrong_j,"j = ",j," - orig_j = ",orig_j
        orig_j=orig_j+n_wrong_j # SKIP WRONG PARAMETERS
        #if orig_j>=57:  # SKIP WRONG PARAMETER
        #  orig_j=orig_j+1       
        #

        reordered_cov[i][j] = orig_cov[orig_i][orig_j]        
        
    outfile.cd()
    reordered_cov   .Write(fit+"_cov")
    reordered_params.Write(fit+"_params")
    
    reordered_covs_d  [fit] = reordered_cov
    reordered_params_d[fit] = reordered_params
    
    
  # Produce untuned central values
  vec = ROOT.TVectorD(n_banff_xsec_params)
  xsec_output_names = ROOT.TObjArray(output_names)
  while xsec_output_names.GetEntries() > n_banff_xsec_params:
    xsec_output_names.RemoveFirst()
    xsec_output_names.Compress()
  xsec_output_names.Write("param_list_banff_xsec",ROOT.TObject.kSingleKey)
  
  for i in xrange(vec.GetNrows()):
    vec[i] = 1.
  for name, val in abnormal_untuned_central_values.iteritems():
    try:
      ind = list(xsec_output_names).index(name)
      vec[ind] = val
    except ValueError:
      print "Making central values - failed to set value of ",name
  
  vec.Write("banff_xsec_untuned_central_values")
    
  infile.Close()
  outfile.Close()
  
  # Start SCC injection
  
  injected_outfile_name = outfile_name.replace(".root","_plus_scc.root")
  injected_outfile = ROOT.TFile(injected_outfile_name,"RECREATE")
  injected_output_names = ROOT.TObjArray()
  for i,p in enumerate(param_order_list):
    injected_output_names.Add(ROOT.TObjString(input_param_list[p]))
    if input_param_list[p] == "MANFFRES":
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
    

    # NEW

    # Check if order is correct (only post-fit, the last matrix stored in arrays)

    if fit == "postfit" and beam_mode_name == "all":
 
      for i in range(len(index_new_orig)):

        new_i  = index_new_orig[i][0]
        orig_i = index_new_orig[i][1]
    
        if orig_i < 0: 
          continue # SCCA and SCCV
        
        if injected_params[new_i] != orig_params[orig_i] :
          print "new index:",new_i,"orig index:",orig_i," --> ",injected_params[new_i],orig_params[orig_i]
          print "New parameters order is wrong!"
          exit(1)

        for j in range(len(index_new_orig)):
            
          new_j  = index_new_orig[j][0]
          orig_j = index_new_orig[j][1]
      
          if orig_j < 0: 
            continue # SCCA and SCCV

          if injected_cov[new_i][new_j] != orig_cov[orig_i][orig_j] :
            print "new index:",new_i,new_j,"orig index:",orig_i,orig_j," --> ",injected_cov[new_i][new_j],orig_cov[orig_i][orig_j]
            print "New matrix order is wrong!"
            exit(1)

      print "\nThe matrix order is OK!\n"


    # Modify covariance elements
    if beam_mode_name=="all" :
      injected_cov[63][63] = 0.0008
      injected_cov[65][65] = 0.0008
      injected_cov[63][65] = -0.0004
      injected_cov[65][63] = -0.0004
    #



    # Write outputs in file

    injected_outfile.cd()
    injected_cov   .Write(fit+"_cov")
    injected_params.Write(fit+"_params")

  # Produce injected untuned central values
  injected_vec = ROOT.TVectorD(n_banff_xsec_params + 2)
  injected_xsec_output_names = ROOT.TObjArray(injected_output_names)
  while injected_xsec_output_names.GetEntries() > n_banff_xsec_params+2:
    injected_xsec_output_names.RemoveFirst()
    injected_xsec_output_names.Compress()
  injected_xsec_output_names.Write("param_list_banff_xsec",ROOT.TObject.kSingleKey)
  
  for i in xrange(injected_vec.GetNrows()):
    injected_vec[i] = 1.
  for name, val in abnormal_untuned_central_values.iteritems():
    try:
      ind = list(injected_xsec_output_names).index(name)
      injected_vec[ind] = val
    except ValueError:
      print "Making SCC injected central values - failed to set value of ",name
  
  injected_vec.Write("banff_xsec_untuned_central_values")
  
  injected_outfile.Close()
