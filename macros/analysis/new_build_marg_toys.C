
//
// ------------------------------------------------------------------------------------------------           
//
// Run toys for marginalization and store into binary files. 
// Toys for different sub-ranges of the free paramters can be run in parallel jobs.
// A file .bin contains the array, while a twin file .bin contains the sub-range of the free parameter
// Each sub-range contains the same number of toys, coming from the same random seed.
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;

TString read_config_file(const char* config_file)
{
   TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
   FileStat_t dummy;

   // Check for config file, and XML data file.
   // If not found in inputs or config, see if it works as a normal path
   if(gSystem->GetPathInfo(config.Data(), dummy)) 
   {
      config = config_file;
      if(gSystem->GetPathInfo(config.Data(), dummy))
      {
         cerr << "Config file not found: " << config_file << endl;
         exit(-1);
      }
   }
   
   return config;
}

void set_global_options(const char* hierarchy, const char* dmsqscheme,
      const char* mixparamchoice, const char* mixparamoctant)
{
   Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
   Opt::Instance()->gOscDmsqScheme = util::DmsqSchemeFromString(dmsqscheme);
   Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
   Opt::Instance()->gOscOctantChoice = util::OctantFromString(mixparamoctant);
   if(Opt::Instance()->Undefined())
   {
      cerr << "Global options not defined correctly. Exiting..." << endl;
      cerr << *Opt::Instance() << endl;
      exit(1);
   }
   Opt::Instance()->SetDefaultParameterValues();
}

void set_mixing_param_choice(double& mixing12_nh, double& mixing13_nh, double& mixing23_nh,
      double& mixing12_ih, double& mixing13_ih, double& mixing23_ih)
{
   bool choose2 = Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta;
   if(mixing13_nh == -99)
   {
      mixing13_nh = choose2 ? conv::kGlobFitSinsq2Theta13 : conv::kGlobFitSinsqTheta13;
   }
   if(mixing23_nh == -99)
   {
      mixing23_nh = choose2 ? conv::kGlobFitSinsq2Theta23 : conv::kGlobFitSinsqTheta23;
   }
   if(mixing12_nh == -99)
   {
      mixing12_nh = choose2 ? conv::kGlobFitSinsq2Theta12 : conv::kGlobFitSinsqTheta12;
   }

   if(mixing13_ih == -99)
   {
      mixing13_ih = choose2 ? conv::kGlobFitSinsq2Theta13 : conv::kGlobFitSinsqTheta13;
   }
   if(mixing23_ih == -99)
   {
      mixing23_ih = choose2 ? conv::kGlobFitSinsq2Theta23 : conv::kGlobFitSinsqTheta23;
   }
   if(mixing12_ih == -99)
   {
      mixing12_ih = choose2 ? conv::kGlobFitSinsq2Theta12 : conv::kGlobFitSinsqTheta12;
   }
}

void set_non_mixing_param(double& dm12sq_nh, double& dm23sq_nh, double& deltacp_nh,
      double& dm12sq_ih, double& dm23sq_ih, double& deltacp_ih, double& beta_bar)
{
   if(dm23sq_nh == -99)
   {
      dm23sq_nh = conv::kGlobFitDm23sq;
   }
   if(dm12sq_nh == -99)
   {
      dm12sq_nh = conv::kGlobFitDm12sq;
   }
   if(deltacp_nh == -99)
   {
      deltacp_nh = conv::kGlobFitDeltaCP;
   }

   if(dm23sq_ih == -99)
   {
      dm23sq_ih = conv::kGlobFitDm23sq;
   }
   if(dm12sq_ih == -99)
   {
      dm12sq_ih = conv::kGlobFitDm12sq;
   }
   if(deltacp_ih == -99)
   {
      deltacp_ih = conv::kGlobFitDeltaCP;
   }

   if(beta_bar == -99)
   {
      beta_bar = conv::kGlobFitBetaBar;
   }
}

PhysParamId_t define_free_parameter(const char* freepar_name)
{
   PhysParamId_t freepar;
   string name(freepar_name);

  // NuBar mode
  if(name == "mixing23Bar")
  {
     freepar = conv::kPhypMixing23Bar;
  }
  else if(name == "mixing13Bar")
  {
     freepar = conv::kPhypMixing13Bar;
  }
  else if(name == "dmsq23Bar")
  {
     freepar = conv::kPhypDm23sqBar;
  }
  else if(name == "dmsq12Bar")
  {
     freepar = conv::kPhypDm12sqBar;
  }
  else if(name == "mixing12Bar")
  {
     freepar = conv::kPhypMixing12Bar;
  }
  else if(name == "deltacpBar")
  {
     freepar= conv::kPhypCPVBar;
  }
  else if(name == "betaBar")
  {
     freepar = conv::kPhypBetaBar;
  }
  // Nu mode
  else if(name =="mixing23")
  {
     freepar = conv::kPhypMixing23;
  }
  else if(name =="mixing13")
  {
     freepar = conv::kPhypMixing13;
  }
  else if(name =="dmsq23")
  {
     freepar = conv::kPhypDm23sq;
  }
  else if(name =="dmsq12")
  {
     freepar = conv::kPhypDm12sq;
  }
  else if(name =="mixing12")
  {
     freepar = conv::kPhypMixing12;
  }
  else if(name =="deltacp")
  {
     freepar = conv::kPhypCPV;
  }
  else
  {
     freepar = conv::kPhypUndefined;
  }

  return freepar;
}

void initialize_parameters(FitParamList* fparm,
      double mixing12, double mixing13, double mixing23,
      double dm12sq, double dm23sq, double deltacp,
      double beta_bar)
{
   conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
   if(hyp == kHypOsc3F)
   {
      fparm->SetInitValue(kPhypMixing12, mixing12);
      fparm->SetInitValue(kPhypMixing13, mixing13);
      fparm->SetInitValue(kPhypMixing23, mixing23);
      fparm->SetInitValue(kPhypCPV, deltacp);
      fparm->SetInitValue(kPhypDm12sq, dm12sq);
      fparm->SetInitValue(kPhypDm23sq, dm23sq);
   }
   else if(hyp == kHypOsc3FBar)
   {
      fparm->SetInitValue(kPhypMixing12, mixing12);
      fparm->SetInitValue(kPhypMixing13, mixing13);
      fparm->SetInitValue(kPhypMixing23, mixing23);
      fparm->SetInitValue(kPhypCPV    ,  deltacp);
      fparm->SetInitValue(kPhypDm12sq,   dm12sq);
      fparm->SetInitValue(kPhypDm23sq,   dm23sq);    
      fparm->SetInitValue(kPhypMixing12Bar, mixing12);
      fparm->SetInitValue(kPhypMixing13Bar, mixing13);
      fparm->SetInitValue(kPhypMixing23Bar, mixing23);
      fparm->SetInitValue(kPhypCPVBar    ,  deltacp);
      fparm->SetInitValue(kPhypDm12sqBar,   dm12sq);
      fparm->SetInitValue(kPhypDm23sqBar,   dm23sq);
   } 
   else if(hyp == kHypOsc3FNormBar)
   {
      fparm->SetInitValue(kPhypMixing12, mixing12);
      fparm->SetInitValue(kPhypMixing13, mixing13);
      fparm->SetInitValue(kPhypMixing23, mixing23);
      fparm->SetInitValue(kPhypCPV    ,  deltacp);
      fparm->SetInitValue(kPhypDm12sq,   dm12sq);
      fparm->SetInitValue(kPhypDm23sq,   dm23sq);    
      fparm->SetInitValue(kPhypBetaBar,  beta_bar);
   }
   else if(hyp == kHypOsc3FBarNormBar) {
      fparm->SetInitValue(kPhypMixing12, mixing12);
      fparm->SetInitValue(kPhypMixing13, mixing13);
      fparm->SetInitValue(kPhypMixing23, mixing23);
      fparm->SetInitValue(kPhypCPV    ,  deltacp);
      fparm->SetInitValue(kPhypDm12sq,   dm12sq);
      fparm->SetInitValue(kPhypDm23sq,   dm23sq);    
      fparm->SetInitValue(kPhypMixing12Bar, mixing12);
      fparm->SetInitValue(kPhypMixing13Bar, mixing13);
      fparm->SetInitValue(kPhypMixing23Bar, mixing23);
      fparm->SetInitValue(kPhypCPVBar    ,  deltacp);
      fparm->SetInitValue(kPhypDm12sqBar,   dm12sq);
      fparm->SetInitValue(kPhypDm23sqBar,   dm23sq);
      fparm->SetInitValue(kPhypBetaBar,     beta_bar);
   }
   else if(hyp == kHypOsc3FBarDissap)
   {
      fparm->SetInitValue(kPhypMixing12, mixing12);
      fparm->SetInitValue(kPhypMixing13, mixing13);
      fparm->SetInitValue(kPhypMixing23, mixing23);
      fparm->SetInitValue(kPhypCPV    ,  deltacp);
      fparm->SetInitValue(kPhypDm12sq,   dm12sq);
      fparm->SetInitValue(kPhypDm23sq,   dm23sq);
      fparm->SetInitValue(kPhypMixing23Bar, mixing23);
      fparm->SetInitValue(kPhypDm23sqBar,   dm23sq);
   }
}

void print_parameters(FitMarginalization* fitmarg)
{
   PhysHypothesis_t hyp = fitmarg->GetFitParamList()->PhysicsHypothesis();
   cout << "   - True oscillation parameters (Normal Hierarchy): " << endl;
   cout << "       > " << util::AsString(kPhypMixing12) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing12) << endl;
   cout << "       > " << util::AsString(kPhypMixing13) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing13) << endl;
   cout << "       > " << util::AsString(kPhypMixing23) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing23) << endl;
   cout << "       > " << util::AsString(kPhypCPV     ) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypCPV     ) << endl;
   cout << "       > " << util::AsString(kPhypDm12sq  ) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm12sq  ) << " eV^2/c^4" << endl;
   cout << "       > " << util::AsString(kPhypDm23sq  ) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm23sq  ) << " eV^2/c^4" << endl;
   if(hyp == kHypOsc3FBar)
   {
      cout << "       > " << util::AsString(kPhypMixing12Bar) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing12Bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing13Bar) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing13Bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing23Bar) << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing23Bar) << endl;
      cout << "       > " << util::AsString(kPhypCPVBar    )  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypCPVBar     ) << endl;
      cout << "       > " << util::AsString(kPhypDm12sqBar)   << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
      cout << "       > " << util::AsString(kPhypDm23sqBar)   << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
   }
   else if(hyp == kHypOsc3FNormBar)
   {
      cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypBetaBar    ) << endl;
   }  
   else if(hyp == kHypOsc3FBarNormBar)
   {
      cout << "       > " << util::AsString(kPhypMixing12Bar)  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing12Bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing13Bar)  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing13Bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing23Bar)  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing23Bar) << endl;
      cout << "       > " << util::AsString(kPhypCPVBar    )   << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypCPVBar     ) << endl;
      cout << "       > " << util::AsString(kPhypDm12sqBar)    << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
      cout << "       > " << util::AsString(kPhypDm23sqBar)    << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
      cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypBetaBar    ) << endl;
   }
   else if(hyp == kHypOsc3FBarDissap)
   {
      cout << "       > " << util::AsString(kPhypMixing23Bar)  << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypMixing23Bar) << endl;
      cout << "       > " << util::AsString(kPhypDm23sqBar)    << " = " << fitmarg->GetFitParamList_NH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
   }

   cout << "   - True oscillation parameters (Inverted Hierarchy): " << endl;
   cout << "       > " << util::AsString(kPhypMixing12) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing12) << endl;
   cout << "       > " << util::AsString(kPhypMixing13) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing13) << endl;
   cout << "       > " << util::AsString(kPhypMixing23) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing23) << endl;
   cout << "       > " << util::AsString(kPhypCPV     ) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypCPV     ) << endl;
   cout << "       > " << util::AsString(kPhypDm12sq  ) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm12sq  ) << " eV^2/c^4" << endl;
   cout << "       > " << util::AsString(kPhypDm23sq  ) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm23sq  ) << " eV^2/c^4" << endl;
   if(hyp == kHypOsc3FBar)
   {
      cout << "       > " << util::AsString(kPhypMixing12Bar) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing12Bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing13Bar) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing13Bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing23Bar) << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing23Bar) << endl;
      cout << "       > " << util::AsString(kPhypCPVBar    )  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypCPVBar     ) << endl;
      cout << "       > " << util::AsString(kPhypDm12sqBar)   << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
      cout << "       > " << util::AsString(kPhypDm23sqBar)   << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
   }
   else if(hyp == kHypOsc3FNormBar)
   {
      cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypBetaBar    ) << endl;
   }  
   else if(hyp == kHypOsc3FBarNormBar)
   {
      cout << "       > " << util::AsString(kPhypMixing12Bar)  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing12Bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing13Bar)  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing13Bar) << endl;
      cout << "       > " << util::AsString(kPhypMixing23Bar)  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing23Bar) << endl;
      cout << "       > " << util::AsString(kPhypCPVBar    )   << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypCPVBar     ) << endl;
      cout << "       > " << util::AsString(kPhypDm12sqBar)    << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm12sqBar  ) << " eV^2/c^4" << endl;
      cout << "       > " << util::AsString(kPhypDm23sqBar)    << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
      cout << "       > " << util::AsString(kPhypBetaBar)      << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypBetaBar    ) << endl;
   }
   else if(hyp == kHypOsc3FBarDissap)
   {
      cout << "       > " << util::AsString(kPhypMixing23Bar)  << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypMixing23Bar) << endl;
      cout << "       > " << util::AsString(kPhypDm23sqBar)    << " = " << fitmarg->GetFitParamList_IH()->InitValue(kPhypDm23sqBar  ) << " eV^2/c^4" << endl;
   }
  
   cout << *fitmarg->GetFitParamList_NH() << endl;
}

void setup_priors(TString file_priors, FitMarginalization* fitmarg, bool reactors)
{
   // Take prior distributions
 
   TString name_prior_th13_dcp_nh = "hPrior_th13_dcp_nh";
   TString name_prior_th13_dcp_ih = "hPrior_th13_dcp_ih";
   TString name_prior_23_nh = "hPrior_L23_nh";
   TString name_prior_23_ih = "hPrior_L23_ih";
 
   TH2D *hPrior_th13_dcp_nh = 0;
   TH2D *hPrior_23_nh = 0;
   TH2D *hPrior_th13_dcp_ih = 0;
   TH2D *hPrior_23_ih = 0;
 
   bool use_prior_dcp = false;
   bool use_prior_23 = false;

   PhysHypothesis_t hyp = fitmarg->GetFitParamList()->PhysicsHypothesis();

   if(file_priors!="")
   {
      use_prior_dcp = true;
      use_prior_23 = true;
    
      cout << endl;
      cout << "Take priors from: " << file_priors << endl;
      TFile *fPrior = new TFile(file_priors.Data(),"READ");
      if(!fPrior->IsOpen() || !fPrior)
      {
         cerr << endl;
         cerr << "Error: the file " << file_priors << " is not open!" << endl;
         cerr << endl;
         exit(1);
      }

      hPrior_th13_dcp_nh = (TH2D*) fPrior->Get(name_prior_th13_dcp_nh.Data());        
      if(hPrior_th13_dcp_nh->IsZombie() || !hPrior_th13_dcp_nh)
      {
         cerr << endl;
         cerr << "Error: the histogram theta13 Vs deltaCP normal hierarchy is not taken!" << endl;
         cerr << endl;
         exit(1);
      }
      hPrior_th13_dcp_ih = (TH2D*) fPrior->Get(name_prior_th13_dcp_ih.Data());       
      if(hPrior_th13_dcp_ih->IsZombie() || !hPrior_th13_dcp_ih)
      {
         cerr << endl;
         cerr << "Error: the histogram theta13 Vs deltaCP inverted hierarchy is not taken!" << endl;
         cerr << endl;
         exit(1);
       }

      hPrior_23_nh = (TH2D*) fPrior->Get(name_prior_23_nh.Data()); 
      if(hPrior_23_nh->IsZombie() || !hPrior_23_nh)
      {
         cerr << endl;
         cerr << "Error: the histogram par23_hist_normal is not taken!" << endl;
         cerr << endl;
         exit(1);
      }
      hPrior_23_ih = (TH2D*) fPrior->Get(name_prior_23_ih.Data()); 
      if(hPrior_23_ih->IsZombie() || !hPrior_23_ih)
      {
         cerr << endl;
         cerr << "Error: the histogram par23_hist_inverted is not taken!" << endl;
         cerr << endl;
         exit(1);
      }
   }
 
   // Set theta13 and deltaCP prior  
   if(use_prior_dcp)
   {
      fitmarg->SetPrior_NH(kPhypCPV,"user",hPrior_th13_dcp_nh);  // mix13 and dcp must have the same prior
      fitmarg->SetPrior_IH(kPhypCPV,"user",hPrior_th13_dcp_ih);

      fitmarg->SetPrior_NH(kPhypMixing13,"user",hPrior_th13_dcp_nh); 
      fitmarg->SetPrior_IH(kPhypMixing13,"user",hPrior_th13_dcp_ih);

      if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar)
      {
         fitmarg->SetPrior_NH(kPhypCPVBar,"user",hPrior_th13_dcp_nh);
         fitmarg->SetPrior_IH(kPhypCPVBar,"user",hPrior_th13_dcp_ih);

         fitmarg->SetPrior_NH(kPhypMixing13Bar,"user",hPrior_th13_dcp_nh);
         fitmarg->SetPrior_IH(kPhypMixing13Bar,"user",hPrior_th13_dcp_ih);
      }
   }
   else
   {  
      // deltaCP
      fitmarg->SetPrior_NH(kPhypCPV, "flat", conv::kPrior_deltaCP_Flat_Low,
		   conv::kPrior_deltaCP_Flat_Up);
      fitmarg->SetPrior_IH(kPhypCPV,"flat", conv::kPrior_deltaCP_Flat_Low,
			conv::kPrior_deltaCP_Flat_Up);
      if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar)
      {
         fitmarg->SetPrior_NH(kPhypCPVBar,"flat", conv::kPrior_deltaCP_Flat_Low,
			   conv::kPrior_deltaCP_Flat_Up);
         fitmarg->SetPrior_IH(kPhypCPVBar,"flat", conv::kPrior_deltaCP_Flat_Low,
			   conv::kPrior_deltaCP_Flat_Up);
      }
      // mixing13
      if(reactors)
      {
         fitmarg->SetPrior_NH(kPhypMixing13,"gaus",
               conv::kPrior_Sinsq2th13_Used_Mean,
      			conv::kPrior_Sinsq2th13_Used_RMS);
         fitmarg->SetPrior_IH(kPhypMixing13,"gaus",
      	      conv::kPrior_Sinsq2th13_Used_Mean,
      		   conv::kPrior_Sinsq2th13_Used_RMS);   
         if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar)
         {
      	   fitmarg->SetPrior_NH(kPhypMixing13Bar,"gaus",
      		   conv::kPrior_Sinsq2th13_Used_Mean,
      			conv::kPrior_Sinsq2th13_Used_RMS);
      	   fitmarg->SetPrior_IH(kPhypMixing13Bar,"gaus",
      			conv::kPrior_Sinsq2th13_Used_Mean,
      			conv::kPrior_Sinsq2th13_Used_RMS);
         }
      }
      else
      {
         fitmarg->SetPrior_NH(kPhypMixing13,"flat",
               conv::kPrior_Sinsq2th13_Flat_Low,
			      conv::kPrior_Sinsq2th13_Flat_Up);
         fitmarg->SetPrior_IH(kPhypMixing13,"flat",
			      conv::kPrior_Sinsq2th13_Flat_Low,
			      conv::kPrior_Sinsq2th13_Flat_Up);   
         if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar)
         {
            fitmarg->SetPrior_NH(kPhypMixing13Bar,"flat",
			      conv::kPrior_Sinsq2th13_Flat_Low,
			      conv::kPrior_Sinsq2th13_Flat_Up);
	         fitmarg->SetPrior_IH(kPhypMixing13Bar,"flat",
			      conv::kPrior_Sinsq2th13_Flat_Low,
			      conv::kPrior_Sinsq2th13_Flat_Up);
         }
      }
   }

   // Set 23 prior
   if(use_prior_23)
   {
      fitmarg->SetPrior_NH(kPhypMixing23,"user",hPrior_23_nh);  // mix23 and dmsq must have the same prior
      fitmarg->SetPrior_IH(kPhypMixing23,"user",hPrior_23_ih);
      fitmarg->SetPrior_NH(kPhypDm23sq,"user",hPrior_23_nh);
      fitmarg->SetPrior_IH(kPhypDm23sq,"user",hPrior_23_ih);
      cout << "setting 23bar priors " << endl;
      if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar || hyp == kHypOsc3FBarDissap)
      {
         cout << "in here" << endl;
         fitmarg->SetPrior_NH(kPhypMixing23Bar,"user",hPrior_23_nh);
         fitmarg->SetPrior_IH(kPhypMixing23Bar,"user",hPrior_23_ih);
         fitmarg->SetPrior_NH(kPhypDm23sqBar,"user",hPrior_23_nh);
         fitmarg->SetPrior_IH(kPhypDm23sqBar,"user",hPrior_23_ih);
      }
   }
   else
   {
      // mixing23 
      fitmarg->SetPrior_NH(kPhypMixing23,"flat",
		   conv::kPrior_Sinsq2th23_Flat_Low,
			conv::kPrior_Sinsq2th23_Flat_Up);
      fitmarg->SetPrior_IH(kPhypMixing23,"flat",
			conv::kPrior_Sinsq2th23_Flat_Low,
			conv::kPrior_Sinsq2th23_Flat_Up);
      cout << "setting 23bar priors " << endl;
      if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar || hyp == kHypOsc3FBarDissap)
      {
         cout << "in here" << endl;
         fitmarg->SetPrior_NH(kPhypMixing23Bar,"flat",
			   conv::kPrior_Sinsq2th23_Flat_Low,
			   conv::kPrior_Sinsq2th23_Flat_Up);
         fitmarg->SetPrior_IH(kPhypMixing23Bar,"flat",
			   conv::kPrior_Sinsq2th23_Flat_Low,
			   conv::kPrior_Sinsq2th23_Flat_Up);
      }
    
      // dmsq23
      fitmarg->SetPrior_NH(kPhypDm23sq,"flat",
            conv::kPrior_Dm23sq_Flat_Low,
			   conv::kPrior_Dm23sq_Flat_Up);
      fitmarg->SetPrior_IH(kPhypDm23sq,"flat",
            conv::kPrior_Dm23sq_Flat_Low,
			   conv::kPrior_Dm23sq_Flat_Up);
      cout << "setting 23bar priors " << endl;
      if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar || hyp == kHypOsc3FBarDissap)
      {
         cout << "in here" << endl;
         fitmarg->SetPrior_NH(kPhypDm23sqBar,"flat",
			   conv::kPrior_Dm23sq_Flat_Low,
			   conv::kPrior_Dm23sq_Flat_Up);
         fitmarg->SetPrior_IH(kPhypDm23sqBar,"flat",
			   conv::kPrior_Dm23sq_Flat_Low,
			   conv::kPrior_Dm23sq_Flat_Up);
      }
   }
  
   // dmsq12
   fitmarg->SetPrior_NH(kPhypDm12sq,"gaus",
         conv::kPrior_Dm12sq_Used_Mean,
		   conv::kPrior_Dm12sq_Used_RMS);
   fitmarg->SetPrior_IH(kPhypDm12sq,"gaus",
         conv::kPrior_Dm12sq_Used_Mean,
		   conv::kPrior_Dm12sq_Used_RMS);

   // mixing12
   fitmarg->SetPrior_NH(kPhypMixing12,"gaus",
		   conv::kPrior_Sinsq2th12_Used_Mean,
		   conv::kPrior_Sinsq2th12_Used_RMS); // prior on sinsq2th
   fitmarg->SetPrior_IH(kPhypMixing12,"gaus",
		   conv::kPrior_Sinsq2th12_Used_Mean,
		   conv::kPrior_Sinsq2th12_Used_RMS);
}

void admin(FitMarginalization* fitmarg, const int fileid, const char* tag)
{
   // Draw prior and toy histograms for marginalized oscillation parameters
   TH1D* hToyHisto_MH;
   TH1D* hPriorHisto_MH;
   PhysHypothesis ph(fitmarg->GetFitParamList()->PhysicsHypothesis());

   const int Nhistos = 8;
   const char* canvas_names_nh[Nhistos] = {
      "cToyHisto_Mix12_nh", "cToyHisto_Dm12_nh", "cToyHisto_Mix13_nh",
      "cToyHisto_dcp_nh", "cToyHisto_Mix23_nh", "cToyHisto_Dm23_nh",
      "cToyHisto_Mix23Bar_nh", "cToyHisto_Dm23Bar_nh"
   };
   const char* canvas_names_ih[Nhistos] = {
      "cToyHisto_Mix12_ih", "cToyHisto_Dm12_ih", "cToyHisto_Mix13_ih",
      "cToyHisto_dcp_ih", "cToyHisto_Mix23_ih", "cToyHisto_Dm23_ih",
      "cToyHisto_Mix23Bar_ih", "cToyHisto_Dm23Bar_ih"
//      "cPriorHisto_MH", "cToyHisto_MH"
   };
   const char* histo_names_nh[Nhistos] = {
      "hToyHisto_Mix12_nh", "hToyHisto_Dm12_nh", "hToyHisto_Mix13_nh",
      "hToyHisto_dcp_nh", "hToyHisto_Mix23_nh", "hToyHisto_Dm23_nh",
      "hToyHisto_Mix23Bar_nh", "hToyHisto_Dm23Bar_nh"
   };
   const char* histo_names_ih[Nhistos] = {
      "hToyHisto_Mix12_ih", "hToyHisto_Dm12_ih", "hToyHisto_Mix13_ih",
      "hToyHisto_dcp_ih", "hToyHisto_Mix23_ih", "hToyHisto_Dm23_ih",
      "hToyHisto_Mix23Bar_ih", "hToyHisto_Dm23Bar_ih"
   };
   const PhysParamId_t params[Nhistos] = {
      conv::kPhypMixing12, conv::kPhypDm12sq, conv::kPhypMixing13, conv::kPhypCPV,
      conv::kPhypMixing23, conv::kPhypDm23sq, conv::kPhypMixing23Bar, conv::kPhypDm23sqBar
   };
   TH1D* toy_histos_nh[Nhistos];
   TCanvas* canvases_nh[Nhistos];

   cout << "Write histograms in output file" << endl;
   TFile* fout = new TFile(Form("%s_histotoys_file%d.root",tag,fileid), "RECREATE");

   // NH
   for(int i = 0; i < Nhistos; ++i)
   {
      toy_histos_nh[i] = 0;
      if(ph.Contains(params[i]))
      {
         TH1D* toy_histo = fitmarg->GetFitParamList_NH()->GetToyHisto(params[i]);
         if(toy_histo)
         {
            canvases_nh[i] = new TCanvas(canvas_names_nh[i], canvas_names_nh[i]);
            toy_histos_nh[i] = new TH1D(*toy_histo);
            toy_histos_nh[i]->Draw("");
            toy_histos_nh[i]->SetName(histo_names_ih[i]);
            toy_histos_nh[i]->Write();
         }
      }
   }

   // IH
   TH1D* toy_histos_ih[Nhistos];
   TCanvas* canvases_ih[Nhistos];

   for(int i = 0; i < Nhistos; ++i)
   {
      toy_histos_ih[i] = 0;
      if(ph.Contains(params[i]))
      {
         TH1D* toy_histo = fitmarg->GetFitParamList_IH()->GetToyHisto(params[i]);
         if(toy_histo)
         {
            canvases_ih[i] = new TCanvas(canvas_names_ih[i], canvas_names_ih[i]);
            toy_histos_ih[i] = new TH1D(*toy_histo);
            toy_histos_ih[i]->Draw("");
            toy_histos_ih[i]->SetName(histo_names_ih[i]);
            toy_histos_ih[i]->Write();
         }
      }
   }

   if(fitmarg->GetHistoPriorMH())
   {
      hPriorHisto_MH = new TH1D(*fitmarg->GetHistoPriorMH());
      TCanvas *cPriorHisto_MH = new TCanvas("cPriorHisto_MH","cPriorHisto_MH");
      hPriorHisto_MH->Draw();
      hPriorHisto_MH->SetName("hPriorHisto_MH");
      hPriorHisto_MH->Write();
   }

   if(fitmarg->GetHistoToyMH())
   {
      hToyHisto_MH = new TH1D(*fitmarg->GetHistoToyMH());
      TCanvas *cToyHisto_MH = new TCanvas("cToyHisto_MH","cToyHisto_MH");
      hToyHisto_MH->Draw();
      hToyHisto_MH->SetName("hToyHisto_MH");
      hToyHisto_MH->Write();
   }

   fitmarg->SaveTreeOfToyParamLists(fout);
   fout->Close();
}

void new_build_marg_toys
(
   const int fileid = 0,
   const int first_bin_freepar0 = -999,
   const int last_bin_freepar0 = -999,
   const int first_bin_freepar1 = -999,
   const int last_bin_freepar1 = -999,
   const char* tag = "PROVA_SystAll_TemplVArray",
   const char* config_file = "analysis_configuration.xml",
   const char* fit_param_list = "nuosc3f+nuisance",
   double input_mixing12_nh = 0.304,
   double input_mixing13_nh = 0.0243,
   double input_mixing23_nh = 0.528,
   double input_deltacp_nh = -1.601,
   double input_dm12sq_nh = 7.53E-5,
   double input_dm23sq_nh = 2.509E-3,
   double input_mixing12_ih = 0.304,
   double input_mixing13_ih = 0.0243,
   double input_mixing23_ih = 0.528,
   double input_deltacp_ih = -1.601,
   double input_dm12sq_ih = 7.53E-5,
   double input_dm23sq_ih = 2.509E-3,
   double input_beta_bar = -99,
   const char* freeparname0 = "mixing13",
   const double freeparmin0 = 0.,
   const double freeparmax0 = 0.4,
   const int NptScan0 = 2,
   const char* freeparname1 = "deltacp",
   const double freeparmin1 = -3.14159265358979312e+00,
   const double freeparmax1 = +3.14159265358979312e+00,
   const int NptScan1 = 2,
   const int NtoysMarg = 100,
   const int tweak_systematics = 1,
   const bool doMargMH = 0,
   const double PriorProbNH = 0.5,
   const bool doRethrowUnphysical = true,
   const char* dmsqscheme = "auto",             // auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi
   const char* hierarchy = "normal",            // normal or inverted
   const char* mixparamchoice = "sinsqtheta",   // sinsq2theta or sinsqtheta
   const char* mixparamoctant = "two",          // one==kNuOctantOne two==kNuOctantTwo
   int seed = 999,
   bool reactors = false
)
{  
   TString file_priors = "";
   TString config = read_config_file(config_file);

   set_global_options(hierarchy, dmsqscheme, mixparamchoice, mixparamoctant);

   // Get superk spectrum predictor & configure
   SuperkExpectation * skexp = new SuperkExpectation;
   skexp->Read(config.Data());

   PhysParamId_t freepar0 = define_free_parameter(freeparname0);
   if(freepar0 == conv::kPhypUndefined)
   {  // Parameter 0 must be always defined
      cout << "Parameter 0 is not defined" << endl;
      exit(1);
   }

   PhysParamId_t freepar1 = define_free_parameter(freeparname1);
   if(freepar1 == conv::kPhypUndefined)
   {
      cout << "Parameter 1 is not defined" << endl;
      cout << "1D fit" << endl;
   }

   // Implement the constructor giving the following input parameters:
   // free param, ntoys, # pt to scan, min par, max par.
   // If min and max not passed, use default ones for that parameter (see FitMarginalization.h)
  
   // Set free param, ntoys, # pt to scan, min par, max par (if min and max not passed, use default ones for that par)

   // free parameter 
   PhysParamId_t parfree[2];
   parfree[0] = freepar0;
   parfree[1] = freepar1;
  
   // # of toys for marginalization
   int ntoys = NtoysMarg;

   // # of free parameter pts to scan 
   int nbinspar[2];
   nbinspar[0] = NptScan0;
   nbinspar[1] = NptScan1;

   // free parameter range
   double minpar[2];
   double maxpar[2];
   minpar[0] = freeparmin0;
   maxpar[0] = freeparmax0;
   minpar[1] = freeparmin1;
   maxpar[1] = freeparmax1;

   // free parameter sub-range 
   int subminbin[2];
   int submaxbin[2];
   subminbin[0] = first_bin_freepar0; submaxbin[0] = last_bin_freepar0;
   subminbin[1] = first_bin_freepar1; submaxbin[1] = last_bin_freepar1;

   // Implement the constructor giving the following input parameters:
   FitMarginalization *fitmarg = new FitMarginalization();   
   fitmarg->SetPhysParFree(parfree,ntoys,nbinspar,minpar,maxpar);
   fitmarg->SetSubRange(subminbin,submaxbin); // By default it scans all the range
   fitmarg->Init(skexp); // Before need SetPhysParFree() to initialize the arrays
   fitmarg->SetSeed(seed);

   // Set fit parameter list
   fitmarg->SetFitParamList(fit_param_list);
   fitmarg->PrintPhysHypothesis();
   //Saving param lists test
   fitmarg->SetSaveParamLists(1);
   //Set marginalization oscillation parameters to defaults  

   FitParamList *fparm    = fitmarg->GetFitParamList();
   FitParamList *fparm_nh = fitmarg->GetFitParamList_NH();
   FitParamList *fparm_ih = fitmarg->GetFitParamList_IH();
   PhysHypothesis_t hyp   = fparm->PhysicsHypothesis();

   set_mixing_param_choice(input_mixing12_nh, input_mixing13_nh, input_mixing23_nh,
         input_mixing12_ih, input_mixing13_ih, input_mixing23_ih);

   set_non_mixing_param(input_dm12sq_nh, input_dm23sq_nh, input_deltacp_nh,
         input_dm12sq_ih, input_dm23sq_ih, input_deltacp_ih, input_beta_bar);
   // NH init param
   initialize_parameters(fparm_nh, input_mixing12_nh, input_mixing13_nh,
         input_mixing23_nh, input_dm12sq_nh, input_dm23sq_nh, input_deltacp_nh,
         input_beta_bar);
   // IH init param
   initialize_parameters(fparm_ih, input_mixing12_ih, input_mixing13_ih,
         input_mixing23_ih, input_dm12sq_ih, input_dm23sq_ih, input_deltacp_ih,
         input_beta_bar);

   if(fparm->PhysicsHypothesis() == kHypNuMuDisapp2F)
   {
      fparm->SetInitValue(kPhypMixing23, input_mixing23);
      fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
   }
 
   print_parameters(fitmarg);

   // Set Prior for Mass Hierarchy (set it for NH and it's automatically done for IH as well)
   fitmarg->SetProbNormalHierarchy(PriorProbNH);
  
   fitmarg->SetMargMassHierarchy(doMargMH);
  
   fitmarg->SetRethrowSystUnphysical(doRethrowUnphysical);

   // Initialize the tree to store the marg toy values of oscillation parameters
   fitmarg->InitTreeOscToys();

   setup_priors(file_priors, fitmarg, reactors);

   // Make toys and store in binary files
    
   // Produce toys for different sub-ranges of the free parameter and store in different files.
   // Useful for running on cluster when many toys are needed (e.g. 10k x 200)
   fitmarg->MakeToysMarginalization(tweak_systematics,skexp); 
   fitmarg->WriteToys2Binary(fileid,tag);  

   // Write tree with stored oscillation parameter values 
   fitmarg->WriteOscToys2Tree(fileid,tag);
  
   admin(fitmarg, fileid, tag);

   return;
}
