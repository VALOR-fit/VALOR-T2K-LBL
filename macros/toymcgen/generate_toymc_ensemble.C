 //
// ------------------------------------------------------------------------------------------------
// Generate / store toy-MC experiments.
// The data can be stored in an XML and/or ROOT format.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Build.h"
#include "AnalysisOptions.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "T2KPhysicalParameterLimits.h"

#include "Utilities.h"
#include "Binning.h"
#include "Conventions/PhysHypothesisId.h"

#include "UtlForMacros.h"
#include "OscParams.h"

#include "TRandom3.h"
#include "TSystem.h"
#include "TStopwatch.h"

#include <sstream>
#include <fstream>
#endif

using namespace std;
using namespace valor;
using namespace valor::conv;

#ifndef __USE_2D_TEMPLATES__ //1D templates
void AddExtraEventsToHist(TH1D * h, string * full_points_to_add);
void AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH1D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand);
#else //2D templates
void AddExtraEventsToHist(TH2D * h, string * full_points_to_add);
void AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH2D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand);
#endif
string * ConstructPointsToAddFromXML(TXMLEngine & xml, XMLNodePointer_t & root_node, valor::conv::SuperkSample_t s, string & sextra, int datasetid);

void generate_toymc_ensemble
(
 const char * tag                     = "toymc",                        /* */
 int          run                     = 0,                              /* */
 double       first_experiment        = 0,                
 double       last_experiment         = 10000,
 const char * hypothesis              = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc, 3f_oscbar" */

 //-999 randomize
 //-99 use global best fits
 double       input_mixing12_nh       = -999,   /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13_nh       = -999,   /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23_nh       = -999,   /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp_nh        = -3.14,   /* deltacp or sindeltacp depending on user choice */
 double       input_dm12sq_nh         = -999,   /* Dm^2_{12},  */
 double       input_dm23sq_nh         = -999,   /* |Dm^2_{23}| */
 double       input_mixing12_ih       = -999,   /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing13_ih       = -999,   /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_mixing23_ih       = -999,   /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp_ih        = -3.14,   /* deltacp or sindeltacp depending on user choice */
 double       input_dm12sq_ih         = -999,   /* Dm^2_{12},  */
 double       input_dm23sq_ih         = -999,   /* |Dm^2_{23}| */
 double       input_beta_bar          = 1.,
 
 // double       input_mixing12_nh       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13_nh       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23_nh       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */   
 // double       input_deltacp_nh        = -1.601, //-99,    /* deltacp, */            
 // double       input_dm12sq_nh         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq_nh         = 2.509E-3, //-99,    /* |Dm^2_{23}| */ 
 // double       input_mixing12_ih       = 0.304,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing13_ih       = 0.0243,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice */
 // double       input_mixing23_ih       = 0.528, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */   
 // double       input_deltacp_ih        = -1.601, //-99,    /* deltacp, */            
 // double       input_dm12sq_ih         = 7.53E-5,    /* Dm^2_{12},  */
 // double       input_dm23sq_ih         = 2.509E-3, //-99,    /* |Dm^2_{23}| */ 
 // double       input_beta_bar          = 1.,

 double       input_prior_prob_nh     = 0.5,   /* prior probability of normal mass hierarchy */
 const bool   doMargMH                = false,  // marginalize Mass Hierarchy?
 
 //string       file_priors             = "priors/joint_T2Kreact_run14/Priors_T2Kreact_run14_Mach3_22Maj15.root",    // file to get priors from (in $T2KOSCANA/inputs) Set to "" for no priors
 string       file_priors             = "",    // file to get priors from (in $T2KOSCANA/inputs) Set to "" for no priors
 
 
 const char * config_file             = "analysis_configuration_oa2017_run1-8_prefit.xml",

 bool rethrow_if_parameters_unphysical = true, // rethrow unphysical toys

 bool         stat_fluctuations       = false,                           /* include statistical fluctuation in toy MC expts? */
 int          tweak_systematics       = false,                              /* generate toy MC with randomized systematics? Numbers greater than 1 are special options */

 string       xml_file_to_add         = "",           /* xml with data in for when creating hybrid data/toy-MC experiments, eg data.xml */
 const int    xml_to_add_dataset_id   = 0,            /* the dataset ID of the xml file to add values from */

 bool         save_fluct_in_xml       = true,                          /* save the values of the variations in the xml file? */

 const char * dmsqscheme              = "auto",                        /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy               = "normal",                       /* normal or inverted */
 const char * mixparamchoice          = "sinsqtheta",                   /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant          = "one",                          /* one==kNuOctantOne two==kNuOctantTwo */
 const char * cpvphase                = "delta",                         /* delta==kNuDeltaCPPhase sindelta==kNuDeltaCPSin */
 bool         reset_xsec_values       = false,
 string       reweight_fake_data      = "MEC2p2hPDD_2017"
)
{
  cout << "start" << endl;
  TStopwatch sw;
  const int randomseed = 100;
  TRandom3 rand(randomseed); //random number to see where inside a bin a toy is generated

  string       points_to_add_1Rmu      = "";                              /* comma separated-list of data-points to add (separate multiple variables with ':'), eg when creating hybrid data/toy-MC experiments */
  string       points_to_add_1Re       = "";                              /* comma separated-list of data-points to add (separate multiple variables with ':'), eg when creating hybrid data/toy-MC experiments */
  string       points_to_add_1Rmu_RHC  = "";                              /* comma separated-list of data-points to add (separate multiple variables with ':'), eg when creating hybrid data/toy-MC experiments */
  string       points_to_add_1Re_RHC   = "";                              /* comma separated-list of data-points to add (separate multiple variables with ':'), eg when creating hybrid data/toy-MC experiments */
  string       points_to_add_MultiRmu  = "";                              /* comma separated-list of data-points to add (separate multiple variables with ':'), eg when creating hybrid data/toy-MC experiments */
  string       points_to_add_MultiRe   = "";                              /* comma separated-list of data-points to add (separate multiple variables with ':'), eg when creating hybrid data/toy-MC experiments */

  
  //some options set here because of CINT 40-argument limit
  bool printsystvariation = true;  /* cout the skexp tweaks */
  bool save_in_xml_format = true;  /* also save in txt format */
  bool save_in_root_format = true;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  string sconfig(config);
  
  // Systematic tweak mode (as string)
  string tweak_sys_str = util::SysThrowIdToString(tweak_systematics);

  // Defined # of experiments
  double nexperiments = last_experiment - first_experiment + 1; 

  // Print out job summary
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Input physics hypothesis: " << hypothesis << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl;
  cout << "   - Randomize systematics : "      << tweak_sys_str << endl;
  cout << "   - Rethrow unphysical : " << ((rethrow_if_parameters_unphysical) ? "YES " : "NO") << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ? "YES " : "NO") << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing toy-MC expts in XML format: "   << ((save_in_xml_format) ?   "YES " : "NO") << endl;
  cout << "     - Storing toy-MC expts in ROOT format:  " << ((save_in_root_format) ?  "YES " : "NO") << endl;

  // Global options
  Opt::Instance()->SetValues(dmsqscheme, hierarchy, mixparamchoice, mixparamoctant, cpvphase);

  //get the prior histograms from the file
  bool UsePrior_th13_dcp = false;
  bool UsePrior_23 = false;
  string file_priors_full_path = "";
  if(!file_priors.empty()) {
    UsePrior_th13_dcp = true;
    UsePrior_23 = true;
    file_priors_full_path = util::FormFilePath(file_priors.c_str(), "inputs");
  }

  


  // Set input osc params to defaults, or set them to be randomised
  util::OscParams::Instance()->SetValues(input_mixing12_nh, input_mixing13_nh, input_mixing23_nh,
                                         input_deltacp_nh,  input_dm12sq_nh,   input_dm23sq_nh,
                                         input_beta_bar,
                                         input_mixing12_ih, input_mixing13_ih, input_mixing23_ih,
                                         input_deltacp_ih,  input_dm12sq_ih,   input_dm23sq_ih,
                                         input_prior_prob_nh, doMargMH,
                                         input_mixing12_nh, input_mixing13_nh, input_mixing23_nh,
                                         input_deltacp_nh,  input_dm12sq_nh,   input_dm23sq_nh,
                                         input_mixing12_ih, input_mixing13_ih, input_mixing23_ih,
                                         input_deltacp_ih,  input_dm12sq_ih,   input_dm23sq_ih);

  // Set priors of osc params to defaults, or set them to be randomised
  util::OscParams::Instance()->SetPrior(conv::kPhypMixing23, "flat", 0.3,          0.7,         "sinsqth" );      
  util::OscParams::Instance()->SetPrior(conv::kPhypDm23sq,   "flat", 2E-3,         3E-3,        ""        );  
  util::OscParams::Instance()->SetPrior(conv::kPhypMixing13, "gaus", 0.093,        0.008,       "sinsq2th"); // Convert to sinsqth
  util::OscParams::Instance()->SetPrior(conv::kPhypCPV,      "flat", -TMath::Pi(), +TMath::Pi(),""        ); 
  util::OscParams::Instance()->SetPrior(conv::kPhypMixing12, "gaus", 0.846,        0.021,       "sinsq2th"); // Convert to sinsqth
  util::OscParams::Instance()->SetPrior(conv::kPhypDm12sq,   "gaus", 7.53E-5,      0.18E-5,     ""        ); 
  

  //util::OscParams::Instance()->Print(cout);


  // Extra points to add to the toy dataset(s)
  string * full_points_to_add[kns] = { 0,0,0,0,0,0 };
  // Are we adding events from an XML file?
  if (!xml_file_to_add.empty()) {
    TXMLEngine xml;
    
    const char * xml_path = util::FormFilePath(xml_file_to_add.c_str(), "inputs");
    XMLDocPointer_t xmldoc = xml.ParseFile(xml_path);
    if (xmldoc == 0) {
      cerr << "Can not parse XML file to add: " << xml_path << endl;
      exit(1);
    }
    XMLNodePointer_t root_node = xml.DocGetRootElement(xmldoc);

    full_points_to_add[0] = ConstructPointsToAddFromXML(xml, root_node, kSKsFCFV1RingMuLike,     points_to_add_1Rmu,     xml_to_add_dataset_id);
    full_points_to_add[1] = ConstructPointsToAddFromXML(xml, root_node, kSKsFCFV1RingELike,      points_to_add_1Re,      xml_to_add_dataset_id);
    full_points_to_add[2] = ConstructPointsToAddFromXML(xml, root_node, kSKsRHCFCFV1RingMuLike,  points_to_add_1Rmu_RHC, xml_to_add_dataset_id);
    full_points_to_add[3] = ConstructPointsToAddFromXML(xml, root_node, kSKsRHCFCFV1RingELike,   points_to_add_1Re_RHC,  xml_to_add_dataset_id);
    full_points_to_add[4] = ConstructPointsToAddFromXML(xml, root_node, kSKsFCFVMultiRingMuLike, points_to_add_MultiRmu, xml_to_add_dataset_id);
    full_points_to_add[5] = ConstructPointsToAddFromXML(xml, root_node, kSKsFCFVMultiRingELike,  points_to_add_MultiRe,  xml_to_add_dataset_id);
  }//xml_file_to_add
  else {
    if ( ! points_to_add_1Rmu.empty()     ) full_points_to_add[0] = new string(points_to_add_1Rmu    );
    if ( ! points_to_add_1Re.empty()      ) full_points_to_add[1] = new string(points_to_add_1Re     );
    if ( ! points_to_add_1Rmu_RHC.empty() ) full_points_to_add[2] = new string(points_to_add_1Rmu_RHC);
    if ( ! points_to_add_1Re_RHC.empty()  ) full_points_to_add[3] = new string(points_to_add_1Re_RHC );
    if ( ! points_to_add_MultiRmu.empty() ) full_points_to_add[4] = new string(points_to_add_MultiRmu);
    if ( ! points_to_add_MultiRe.empty()  ) full_points_to_add[5] = new string(points_to_add_MultiRe );
  }//!xml_file_to_add
  cout << "   - Adding the following data-points in each toy-MC expt: " << endl;
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    if (s == conv::kSKsUndefined) continue;
    cout << "    For " << util::AsStringCompact(s) << "\t" << ((full_points_to_add[is]) ? *full_points_to_add[is] : "*None*") << endl;
  }//is
  cout << endl;
  
  // Open files
  TFile * f = 0;
  if(save_in_root_format) {
    ostringstream root_filename;
    root_filename << tag << "-" << run << ".root";
    cout << root_filename.str();
    f = new TFile(root_filename.str().c_str(), "recreate");
  }

  ofstream fxml;
  ofstream ftxt;
  if(save_in_xml_format) {
    ostringstream xml_filename;
    xml_filename << tag << "-" << run << ".xml";
    fxml.open(xml_filename.str().c_str());
    fxml << "<?xml version=\"1.0\"?> " << endl;
    fxml << "<analysis_inputs>       " << endl;
    fxml << "  <superk>              " << endl;
    fxml << "    <observation>       " << endl;
    ostringstream txt_filename;
    txt_filename << tag << "-" << run << ".txt";
    ftxt.open(txt_filename.str().c_str());
  }

  //Making fake data reweight
  if (reweight_fake_data.compare("")) {
    util::UtlMakeFakeData::Instance()->SetReweightStyle(reweight_fake_data);
    util::UtlMakeFakeData::Instance()->Activate();
  }

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(sconfig);
  
  int weird_applied_neg      = 0;
  int weird_unphysical_tweak = 0;
  for(int iexperiment=first_experiment; iexperiment <= last_experiment; iexperiment++) {
     cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;

     // Generate mock data
     cout << "- Generating toy-MC..." << endl;
     SuperkObservation * skobs = new SuperkObservation;
     PhysHypothesis ph(hypothesis);


     // Randomize specified osc params for current toy MC expt
     util::OscParams::Instance()->Randomize();

     // Set the mass hierarchy 
     // Do it after randomizing MH, otherwise we can have two different MH values and it crashes
     util::SetMassHierarchy(skexp); 

     // Set physics hypothesis
     util::OscParams::Instance()->UpdatePhysHypothesis(ph);

     // Fill toy histograms of osc params
     util::OscParams::Instance()->FillOscToyHisto(ph);
     
     //
     
     bool applied_negevents_condition = true;
     bool unphysical_tweak = true;
     while(applied_negevents_condition || unphysical_tweak) {

       util::PerformSysThrow(tweak_systematics, skexp, printsystvariation);
       if(reset_xsec_values)
       {
         //          twkv[pfo_ip] = 1 + ((225. - 225.)/225.) ;
         skexp->GetBANFF()->SetCurrentValue(conv::kNuisBANFFMEC, 1);
         skexp->GetBANFF()->SetCurrentValue(conv::kNuisBANFFMAQE, 1);
         skexp->GetBANFF()->PrintSystematicVariations();
       }
       // Rethrow unphysical parameters
       PhysParamList true_nuisance_params = skexp->Tweaks();
       if (rethrow_if_parameters_unphysical) {
	 unphysical_tweak = !(valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params));
	 if (unphysical_tweak) {
	   cout << "Parameters are unphysical - rethrowing" << std::endl;
	   weird_unphysical_tweak++;
	   continue;
	 }
       }
       else {
	 unphysical_tweak = false;
       }
       if  ( reweight_fake_data =="MEC2p2hPDD_2017") {
         cout << "Setting the value of the 2p2h Shape Parameter to 1 " << endl;
         skexp->GetBANFF()->SetCurrentValue(conv::kNuisBANFFShape2p2h, 1);
       }
       skexp->GetBANFF()->PrintSystematicVariations();
       applied_negevents_condition = skobs->CreateMockData(*skexp, ph, stat_fluctuations);
       cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
       if(applied_negevents_condition)
	 weird_applied_neg++;
     }//applied_negevents_condition

     //Get the spectra
#ifndef __USE_2D_TEMPLATES__ //1D templates
     TH1D * mock_data[kns];
#else
     TH2D * mock_data[kns];
#endif
     for (int is = 0 ; is < conv::kns ; is++) {
       SuperkSample_t s = conv::kSKSamples[is];
       if (s == conv::kSKsUndefined) continue;
#ifndef __USE_2D_TEMPLATES__ //1D templates
       mock_data[is] = skobs->ErecoSpectrum(s);
#else
       mock_data[is] = skobs->TwoDSpectrum(s);
#endif
       mock_data[is]->SetDirectory(0);
     }//is
      
     // If an XML file was specified, add an entry
     if(save_in_xml_format) {

       // Dataset info
       fxml << "      " << Form("<data_set_%d>",iexperiment) << endl;
       if(save_fluct_in_xml) {
         fxml << "<!-- ";
         ph.Print(fxml);
         if(tweak_systematics)
           skexp->PrintSystematicVariations(fxml);
         fxml << " -->" << endl;
       }

       for (int is = 0 ; is < conv::kns ; is++) {
	 SuperkSample_t s = conv::kSKSamples[is];
	 if (s == conv::kSKsUndefined) continue;
	 AddEventsToXML(fxml, ftxt, mock_data[is], s, full_points_to_add[is], rand);
       }//is
       
       fxml << "      " << Form("</data_set_%d>\n",iexperiment) << endl;       
     }//save to xml file

     // If a ROOT filename was specified, save the ROOT histogram too
     if(f) {
       for (int is = 0 ; is < conv::kns ; is++) {
	 SuperkSample_t s = conv::kSKSamples[is];
	 if (s == conv::kSKsUndefined) continue;
	 AddExtraEventsToHist(mock_data[is], full_points_to_add[is]);
       }//is
      
       f->cd();

       for (int is = 0 ; is < conv::kns ; is++) {
	 SuperkSample_t s = conv::kSKSamples[is];
	 if (s == conv::kSKsUndefined) continue;
	 mock_data[is]->Write(Form("mock_data_%s_%d", util::AsStringCompact(s).c_str(), iexperiment));
	 delete mock_data[is];
       }//is
            
       // Object for storing input oscillation and nuisance params in each toy-MC
       PhysParamList * iparm = new PhysParamList;
       iparm->Clear();
       // keep track of oscillation parameters
       PhysParamList true_oscillation_params = ph.AsPhysParamList();
       iparm->Append(true_oscillation_params);
       // keep track of all systematic params
       PhysParamList true_nuisance_params = skexp->Tweaks();
       iparm->Append(true_nuisance_params);
       // save it
       iparm->Write(Form("iparm_%d",iexperiment));
       delete iparm;
     }//save to ROOT file

  }//iexperiment


  // Store toy histograms    
  char name[200];
  sprintf(name,"%s_oscparam_toys.root",tag);
  TFile *ftoys = new TFile(name,"recreate");
  util::OscParams::Instance()->Write(ftoys);
  ftoys->Close();

  //
  // Close files
  //
  if(save_in_xml_format) {
    fxml << "    </observation>       " << endl;
    fxml << "  </superk>              " << endl;
    fxml << "</analysis_inputs>       ";// << endl;
    
    fxml << endl;
    fxml << "<!--" << endl;
    util::OscParams::Instance()->Print_Toys(fxml);
    fxml << "-->" << endl;
    fxml << endl;

    fxml.close();
  }
  if(f) {
    f->Close();
    delete f;
  }//f
  
  cout << endl 
       << " - Summary - " << endl
       << "    * Created " << nexperiments << " toys" << endl
       << "    * " << weird_unphysical_tweak << " unphysical_tweaks rethrown" << endl
       << "    * " << weird_applied_neg << " applied_negevents_condition rethrown" << endl
       << endl;

  sw.Stop();
  sw.Print();
}











#ifndef __USE_2D_TEMPLATES__ //1D templates
void AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH1D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand)
#else
void AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH2D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand)
#endif
{
  ftxt << util::AsStringCompact(s) << ":" << endl;
  fxml << "        <" << util::GetXMLDatasetName(s) << ">" << endl;
  fxml << "          <" << util::GetXMLDatasetParamsName() << ">" << endl;

  int ievtot = 0;

#ifndef __USE_2D_TEMPLATES__ //1D templates
  const int n_points_per_line_xml = 10; // this should be an option. But 1D templates are going so there's no point...

  // Adding Ereco values as a comma separated list.
  // Assuming uniform distribution within each bin, randomize values within each bin.
  for(int ibin = 1; ibin <= h->GetNbinsX(); ibin++) {
    double Nb   = h->GetBinContent (ibin); // number of events in this Ereco bin
    double Emin = h->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
    double Emax = h->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
    //double Ecent = h->GetXaxis()->GetBinCenter (ibin); // Ereco bin center
    for(int iev = 0; iev < Nb; iev++) {
      double E = Emin + rand.Uniform() * (Emax - Emin);
      if(ievtot > 0) {
	fxml << ", "; 
	ftxt << ","; 
	if ( ! (ievtot % n_points_per_line_xml) ) fxml << endl; 
      }
      fxml << E; 
      ftxt << E; 
      ievtot++;
    }//iev
  }//ibin
  // add extra data-points, if any
  if(full_points_to_add) {
    fxml << ",\n" << endl; // should really have an XML comment here but TXMLEngine won't let us.
    fxml << *full_points_to_add;
    ftxt << "," << *full_points_to_add;
  }
  fxml << endl;

#else //2D templates
       
  // Adding 2D template values with each event as its own node containing 2 values
  // Assuming uniform distribution within each bin, randomize values within each bin.
  for(int ibinx = 1; ibinx <= h->GetNbinsX(); ibinx++) {
    double Xmin  = h->GetXaxis()->GetBinLowEdge(ibinx); // X bin min
    double Xmax  = h->GetXaxis()->GetBinUpEdge (ibinx); // X bin max
    //double Xcent = h->GetXaxis()->GetBinCenter (ibinx); // X bin center
    for(int ibiny = 1; ibiny <= h->GetNbinsY(); ibiny++) {
      double Nb    = h->GetBinContent (ibinx, ibiny); // number of events in this bin
      double Ymin  = h->GetYaxis()->GetBinLowEdge(ibiny); // Y bin min
      double Ymax  = h->GetYaxis()->GetBinUpEdge (ibiny); // Y bin max
      //double Ycent = h->GetYaxis()->GetBinCenter (ibiny); // Y bin center
      for(int iev = 0; iev < Nb; iev++) {
	double X = Xmin + rand.Uniform() * (Xmax - Xmin);
	double Y = Ymin + rand.Uniform() * (Ymax - Ymin);
	if(ievtot > 0) {
	  ftxt << ": "; 
	}
	fxml << "            <ev" << ievtot << ">"
	     << X << "," << Y
	     << "</ev" << ievtot << ">"
	     << endl;
	ftxt << X << "," << Y; 
	ievtot++;
      }//iev
    }//ibiny
  }//ibinx
  // add extra data-points, if any
  if(full_points_to_add) {
    fxml << endl; // should really have an XML comment here but TXMLEngine won't let us.
    vector<string> vecpairs = util::Tokens(*full_points_to_add,",");
    assert(vecpairs.size()>0);
    for(unsigned int ipair = 0; ipair < vecpairs.size(); ipair++) {
      vector<string> values = util::Tokens(vecpairs[ipair], ":");
      assert(values.size() == 2);
      fxml << "            <ev" << ievtot << ">"
	   << atof(values[0].c_str()) << "," << atof(values[1].c_str())
	   << "</ev" << ievtot << ">"
	   << endl;
      ievtot++;
    }//ipair
    ftxt << ":" << *full_points_to_add;
  }
  // add the nevents to the XML
  fxml << "            <nevents> " << ievtot << "</nevents>" << endl;
#endif //__USE_2D_TEMPLATES__
	 
  ftxt << endl;
  fxml << "          </" << util::GetXMLDatasetParamsName() << "> " << endl;
  fxml << "        </" << util::GetXMLDatasetName(s) << ">    " << endl;       
}


#ifndef __USE_2D_TEMPLATES__ //1D templates
void AddExtraEventsToHist(TH1D * h, string * full_points_to_add)
#else
void AddExtraEventsToHist(TH2D * h, string * full_points_to_add)
#endif
{
  if(full_points_to_add) {
#ifndef __USE_2D_TEMPLATES__ //1D templates
    vector<string> vecp = util::Tokens(*full_points_to_add, ",");
    assert(vecp.size() > 0);
    for(unsigned int ip=0; ip < vecp.size(); ip++) {
      h->Fill(atof(vecp[ip].c_str()));
    }//ip
#else //2D templates
    vector<string> vecpairs = util::Tokens(*full_points_to_add,",");
    assert(vecpairs.size()>0);
    for(unsigned int ipair = 0; ipair < vecpairs.size(); ipair++) {
      vector<string> values = util::Tokens(vecpairs[ipair], ":");
      assert(values.size() == 2);
      h->Fill(atof(values[0].c_str()), atof(values[1].c_str()));
    }//ipair
#endif
  }//full_points_to_add
}

string * ConstructPointsToAddFromXML(TXMLEngine & xml, XMLNodePointer_t & root_node, valor::conv::SuperkSample_t s, string & sextra, int datasetid)
{
  if(util::ArrIdx::Instance()->SuperkSampleIdx(s) == -1)
    return 0;
#ifndef __USE_2D_TEMPLATES__ //1D templates
  string sxml = util::XML2String(&xml, root_node, Form("superk/observation/data_set_%d/%s/ereco_values", datasetid, util::GetXMLMCTreeName(s).c_str()));
  string * sout = new string(sxml);
  if(!sextra.empty()) {
    sout->append(",");
    sout->append(sextra);
  }
#else
  int nevents = util::XML2Int(&xml, root_node, Form("superk/observation/data_set_%d/%s/%s/nevents",
						    datasetid, util::GetXMLMCTreeName(s).c_str(), util::GetXMLDatasetParamsName().c_str()));
  string * sout = new string();
  for(int iev = 0; iev < nevents; iev++) {
    if(iev != 0)
      sout->append(", ");
    string sxml = util::XML2String(&xml, root_node, Form("superk/observation/data_set_%d/%s/%s/ev%d",
							 datasetid, util::GetXMLMCTreeName(s).c_str(), util::GetXMLDatasetParamsName().c_str(), iev));
    std::replace(sxml.begin(), sxml.end(), ',', ':');
    sout->append(sxml);
  }//iev
  if(!sextra.empty()) {
    sout->append(", ");
    sout->append(sextra);
  }
#endif
  return sout;
}
