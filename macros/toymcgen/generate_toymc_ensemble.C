//
// ------------------------------------------------------------------------------------------------
// Generate / store toy-MC experiments.
// The data can be stored in an XML and/or ROOT format.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "T2KPhysicalParameterLimits.h"

#include "Utilities.h"
#include "Binning.h"
#include "Utilities/UtlXML.h"
#include "Conventions/PhysHypothesisId.h"

#include "Build.h"
#include "UtlForMacros.h"

#include "AnalysisOptions.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

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
string * ConstructPointsToAddFromXML(TXMLEngine & xml, XMLNodePointer_t & root_node, valor::conv::SuperkSample_t s, string sextra);

void generate_toymc_ensemble
(
 const char * tag                     = "toymc_numudisap",                        /* */
 int          run                     = 0,                              /* */
 double       first_experiment        = 0,                
 double       last_experiment         = 0,
 const char * hypothesis              = "3f_oscbar",                       /* options: "2f_numu_disap", "3f_osc, 3f_oscbar" */

 double       input_mixing12_nh       = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13_nh       = -99,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23_nh       = -99,   //conv::kT2KReacBestFitSinsqTheta23_NH, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_deltacp_nh        = -99,   //conv::kT2KReacBestFitDeltaCP_NH, //-99,    /* deltacp, */
 double       input_dm12sq_nh         = -99,    /* Dm^2_{12},  */
 double       input_dm23sq_nh         = -99,   //conv::kT2KReacBestFitDmsq23_NH, //-99,    /* |Dm^2_{23}| */
 double       input_mixing12_ih       = -99,    /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13_ih       = -99,    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23_ih       = -99,   //conv::kT2KReacBestFitSinsqTheta23_IH, //-99,    /* 23-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_deltacp_ih        = -99,   //conv::kT2KReacBestFitDeltaCP_IH, //-99,    /* deltacp, if -999 then randomize */
 double       input_dm12sq_ih         = -99,    /* Dm^2_{12},  */
 double       input_dm23sq_ih         = -99,   //conv::kT2KReacBestFitDmsq23_IH, //-99,    /* |Dm^2_{23}|, if -999 then randomize */
 double       input_beta_bar          = -99,
 
 const char * config_file      = "analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml",

 bool rethrow_if_parameters_unphysical = true, // rethrow unphysical toys

 bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
 int          tweak_systematics       = 1,                              /* generate toy MC with randomized systematics? Numbers greater than 1 are special options */
 bool         printsystvariation      = true,                           /* cout the skexp tweaks */

 bool         save_in_xml_format      = true,                           /* also save in txt format */
 bool         save_in_root_format     = true,                           /* */  
 const char * points_to_add_1Rmu      = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_1Re       = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_1Rmu_RHC  = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_1Re_RHC   = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_MultiRmu  = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 const char * points_to_add_MultiRe   = 0,                              /* comma separated-list of data-points to add, eg when creating hybrid data/toy-MC experiments */
 //const char * xml_file_to_add         = "data_newstruct.xml",           /* xml with data in for when creating hybrid data/toy-MC experiments, eg data.xml. Always uses set 0 */
 const char * xml_file_to_add         = 0,           /* xml with data in for when creating hybrid data/toy-MC experiments, eg data.xml. Always uses set 0 */
 bool         save_fluct_in_xml       = true,                          /* save the values of the variations in the xml file? */
 const char * dmsqscheme              = "auto",                        /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy               = "normal",                       /* normal or inverted */
 const char * mixparamchoice          = "sinsqtheta",                   /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant          = "one",                          /* one==kNuOctantOne two==kNuOctantTwo */
 const char * cpvphase                = "delta"                         /* delta==kNuDeltaCPPhase sindelta==kNuDeltaCPSin */
 )
{
  TRandom3 rand(0);

  TH1D * h_0 = new TH1D ("h_0", "h_0", 100, 0, 10);
  TH1D * h_1 = new TH1D ("h_1", "h_1", 100, 0, 10);
  TH1D * h_2 = new TH1D ("h_2","h_2", 100, 0, 10);
  TH1D * h_3 = new TH1D ("h_3","h_3", 100, 0, 10);
  TH1D * h_4 = new TH1D ("h_4","h_4", 100, 0, 10);
  TH2D * h_4_Escale = new TH2D ("h_4_Escale","h_4_Escale", 100, -5, 5, 100, 0, 10);

  // Analysis configuration file
  const char * config = util::FormFilePath(config_file, "config");

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
  util::SetGlobalOptions(dmsqscheme, hierarchy, mixparamchoice, mixparamoctant, cpvphase);

  // Set input osc params to defaults, or set them to be randomised
  bool randomize_mixing12_nh, randomize_mixing13_nh, randomize_mixing23_nh;
  bool randomize_deltacp_nh, randomize_dm12sq_nh, randomize_dm23sq_nh;
  bool randomize_beta_bar;
  bool randomize_mixing12_ih, randomize_mixing13_ih, randomize_mixing23_ih;
  bool randomize_deltacp_ih, randomize_dm12sq_ih, randomize_dm23sq_ih;
  util::SetOscParamOptions(input_mixing12_nh, input_mixing13_nh, input_mixing23_nh,
			   input_deltacp_nh,  input_dm12sq_nh,   input_dm23sq_nh,
			   randomize_mixing12_nh, randomize_mixing13_nh, randomize_mixing23_nh,
			   randomize_deltacp_nh,  randomize_dm12sq_nh,   randomize_dm23sq_nh,
			   input_beta_bar, randomize_beta_bar,
			   input_mixing12_ih, input_mixing13_ih, input_mixing23_ih,
			   input_deltacp_ih,  input_dm12sq_ih,   input_dm23sq_ih,
			   randomize_mixing12_ih, randomize_mixing13_ih, randomize_mixing23_ih,
			   randomize_deltacp_ih,  randomize_dm12sq_ih,   randomize_dm23sq_ih		     
			   );
  
  string * full_points_to_add_1Rmu     = 0;
  string * full_points_to_add_1Re      = 0;
  string * full_points_to_add_1Rmu_RHC = 0;
  string * full_points_to_add_1Re_RHC  = 0;
  string * full_points_to_add_MultiRmu = 0;
  string * full_points_to_add_MultiRe  = 0;
  
  // Are we adding events from an XML file?
  if (xml_file_to_add) {
    cerr << "xml_file_to_add needs testing (and not yet implemented for 2D templates)" << endl;
    exit(-1);

    TXMLEngine xml;
    
    const char * xml_path = util::FormFilePath(xml_file_to_add, gSystem->Getenv("T2KOSCANA"));
    XMLDocPointer_t xmldoc = xml.ParseFile(xml_path);
    if (xmldoc==0) {
      cerr << "Can not parse XML file to add: " << xml_path << endl;
      exit(1);
    }
    XMLNodePointer_t root_node = xml.DocGetRootElement(xmldoc);

    full_points_to_add_1Rmu     = ConstructPointsToAddFromXML(xml, root_node, kSKsFCFV1RingMuLike,     points_to_add_1Rmu);
    full_points_to_add_1Re      = ConstructPointsToAddFromXML(xml, root_node, kSKsFCFV1RingELike,      points_to_add_1Re);
    full_points_to_add_1Rmu_RHC = ConstructPointsToAddFromXML(xml, root_node, kSKsRHCFCFV1RingMuLike,  points_to_add_1Rmu_RHC);
    full_points_to_add_1Re_RHC  = ConstructPointsToAddFromXML(xml, root_node, kSKsRHCFCFV1RingELike,   points_to_add_1Re_RHC);
    full_points_to_add_MultiRmu = ConstructPointsToAddFromXML(xml, root_node, kSKsFCFVMultiRingMuLike, points_to_add_MultiRmu);
    full_points_to_add_MultiRe  = ConstructPointsToAddFromXML(xml, root_node, kSKsFCFVMultiRingELike,  points_to_add_MultiRe);
  }//xml_file_to_add
  else {
    if (points_to_add_1Rmu    ) full_points_to_add_1Rmu     = new string(points_to_add_1Rmu    );
    if (points_to_add_1Re     ) full_points_to_add_1Re      = new string(points_to_add_1Re     );
    if (points_to_add_1Rmu_RHC) full_points_to_add_1Rmu_RHC = new string(points_to_add_1Rmu_RHC);
    if (points_to_add_1Re_RHC ) full_points_to_add_1Re_RHC  = new string(points_to_add_1Re_RHC );
    if (points_to_add_MultiRmu) full_points_to_add_MultiRmu = new string(points_to_add_MultiRmu);
    if (points_to_add_MultiRe ) full_points_to_add_MultiRe  = new string(points_to_add_MultiRe );
  }//!xml_file_to_add
  cout << "   - Adding the following data-points in each toy-MC expt: " << endl;
  cout << "    For 1R mu-like: "     << ((full_points_to_add_1Rmu)     ? *full_points_to_add_1Rmu     : "*None*") << endl;
  cout << "    For 1R e-like: "      << ((full_points_to_add_1Re)      ? *full_points_to_add_1Re      : "*None*") << endl;
  cout << "    For 1R mu-like RHC: " << ((full_points_to_add_1Rmu_RHC) ? *full_points_to_add_1Rmu_RHC : "*None*") << endl;
  cout << "    For 1R e-like RHC: "  << ((full_points_to_add_1Re_RHC)  ? *full_points_to_add_1Re_RHC  : "*None*") << endl;
  cout << "    For MultiR mu-like: " << ((full_points_to_add_MultiRmu) ? *full_points_to_add_MultiRmu : "*None*") << endl;
  cout << "    For MultiR e-like: "  << ((full_points_to_add_MultiRe)  ? *full_points_to_add_MultiRe  : "*None*") << endl;
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

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  int weird = 0;
  for(int iexperiment=first_experiment; iexperiment <= last_experiment; iexperiment++) {
     cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;

     // Generate mock data
     cout << "- Generating toy-MC..." << endl;
     SuperkObservation * skobs = new SuperkObservation;
     PhysHypothesis ph(hypothesis);

     // Randomize specified osc params for current toy MC expt
     util::RandomizeOscParams(input_mixing12_nh, input_mixing13_nh, input_mixing23_nh,
			      input_deltacp_nh,  input_dm12sq_nh,   input_dm23sq_nh,
			      randomize_mixing12_nh, randomize_mixing13_nh, randomize_mixing23_nh,
			      randomize_deltacp_nh,  randomize_dm12sq_nh,   randomize_dm23sq_nh,
			      input_beta_bar, randomize_beta_bar,
			      input_mixing12_ih, input_mixing13_ih, input_mixing23_ih,
			      input_deltacp_ih,  input_dm12sq_ih,   input_dm23sq_ih,
			      randomize_mixing12_ih, randomize_mixing13_ih, randomize_mixing23_ih,
			      randomize_deltacp_ih,  randomize_dm12sq_ih,   randomize_dm23sq_ih		     
			      );

     // Set physics hypothesis
     //TODO use priors to decide which hierarchy to set to
     util::SetOscParams(ph,
			input_mixing12_nh, input_mixing13_nh, input_mixing23_nh,
			input_deltacp_nh,  input_dm12sq_nh,   input_dm23sq_nh,
			input_beta_bar
			);





     bool applied_negevents_condition = true;
     bool unphysical_tweak = true;
     while(applied_negevents_condition || unphysical_tweak) {

       util::PerformSysThrow(tweak_systematics, skexp, printsystvariation);

       // Rethrow unphysical parameters
       PhysParamList true_nuisance_params = skexp->Tweaks();
       if (rethrow_if_parameters_unphysical) {
	 unphysical_tweak = !(valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params));
	 if (unphysical_tweak) {
	   cout << "Parameters are unphysical - rethrowing" << std::endl;
	   continue;
	 }
       }
       else {
	 unphysical_tweak = false;
       }

       applied_negevents_condition = skobs->CreateMockData(*skexp, ph, stat_fluctuations);
       cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
     }//applied_negevents_condition
  
     //Get the spectra
#ifdef __FHC_RUNNING__
#ifndef __USE_2D_TEMPLATES__ //1D templates
#ifndef __NUEONLY__
     TH1D * mock_data_1Rmu = skobs->ErecoSpectrum(kSKsFCFV1RingMuLike);
#endif
#ifndef __NUMUONLY__
     TH1D * mock_data_1Re  = skobs->ErecoSpectrum(kSKsFCFV1RingELike);
#endif
#else //2D templates
#ifndef __NUEONLY__
     TH2D * mock_data_1Rmu = skobs->TwoDSpectrum(kSKsFCFV1RingMuLike);
#endif
#ifndef __NUMUONLY__
     TH2D * mock_data_1Re  = skobs->TwoDSpectrum(kSKsFCFV1RingELike);
#endif
#endif //__USE_2D_TEMPLATES__
#ifndef __NUEONLY__
     mock_data_1Rmu->SetDirectory(0);
#endif
#ifndef __NUMUONLY__
     mock_data_1Re->SetDirectory(0);
#endif
#endif //__FHC_RUNNING__
     
#ifdef __RHC_RUNNING__
#ifndef __USE_2D_TEMPLATES__ //1D templates
#ifndef __NUEONLY__
     TH1D * mock_data_1Rmu_RHC = skobs->ErecoSpectrum(kSKsRHCFCFV1RingMuLike);
#endif
#ifndef __NUMUONLY__
     TH1D * mock_data_1Re_RHC  = skobs->ErecoSpectrum(kSKsRHCFCFV1RingELike);
#endif
#else //2D templates
#ifndef __NUEONLY__
     TH2D * mock_data_1Rmu_RHC = skobs->TwoDSpectrum(kSKsRHCFCFV1RingMuLike);
#endif
#ifndef __NUMUONLY__
     TH2D * mock_data_1Re_RHC  = skobs->TwoDSpectrum(kSKsRHCFCFV1RingELike);
#endif
#endif //__USE_2D_TEMPLATES__
#ifndef __NUEONLY__
     mock_data_1Rmu_RHC->SetDirectory(0);
#endif
#ifndef __NUMUONLY__
     mock_data_1Re_RHC->SetDirectory(0);
#endif
#endif //__RHC_RUNNING__
     
#ifdef __MULTIRING_FHC_RUNNING__
#ifndef __USE_2D_TEMPLATES__ //1D templates
#ifndef __NUEONLY__
     TH1D* mock_data_MultiRmu = skobs->ErecoSpectrum(kSKsFCFVMultiRingMuLike);
#endif
#ifndef __NUMUONLY__
     TH1D* mock_data_MultiRe  = skobs->ErecoSpectrum(kSKsFCFVMultiRingELike);
#endif
#else //2D templates
#ifndef __NUEONLY__
     TH2D* mock_data_MultiRmu = skobs->TwoDSpectrum(kSKsFCFVMultiRingMuLike);
#endif
#ifndef __NUMUONLY__
     TH2D* mock_data_MultiRe  = skobs->TwoDSpectrum(kSKsFCFVMultiRingELike);     
#endif
#endif //__USE_2D_TEMPLATES__
#ifndef __NUEONLY__
     mock_data_MultiRmu->SetDirectory(0);     
#endif
#ifndef __NUMUONLY__
     mock_data_MultiRe->SetDirectory(0);
#endif
#endif //__MULTIRING_FHC_RUNNING__
     
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
     
#ifdef __FHC_RUNNING__
#ifndef __NUEONLY__
       // Now the 1Rmu part
       AddEventsToXML(fxml, ftxt, mock_data_1Rmu, kSKsFCFV1RingMuLike, full_points_to_add_1Rmu, rand);
#endif
#ifndef __NUMUONLY__
       //Now the 1Re part
       AddEventsToXML(fxml, ftxt, mock_data_1Re,  kSKsFCFV1RingELike,  full_points_to_add_1Re, rand);
#endif
#endif //__FHC_RUNNING__
       
#ifdef __RHC_RUNNING__
#ifndef __NUEONLY__
       // Now the 1Rmu RHC part
       AddEventsToXML(fxml, ftxt, mock_data_1Rmu_RHC, kSKsRHCFCFV1RingMuLike, full_points_to_add_1Rmu_RHC, rand);
#endif
#ifndef __NUMUONLY__
       //Now the 1Re RHC part
       AddEventsToXML(fxml, ftxt, mock_data_1Re_RHC,  kSKsRHCFCFV1RingELike,  full_points_to_add_1Re_RHC, rand);
#endif
#endif //__RHC_RUNNING__
       
#ifdef __MULTIRING_FHC_RUNNING__
#ifndef __NUEONLY__
       // Now the MultiRmu part
       AddEventsToXML(fxml, ftxt, mock_data_MultiRmu, kSKsFCFVMultiRingMuLike, full_points_to_add_MultiRmu, rand);
#endif
#ifndef __NUMUONLY__
       //Now the 1Re RHC part
       AddEventsToXML(fxml, ftxt, mock_data_MultiRe,  kSKsFCFVMultiRingELike,  full_points_to_add_MultiRe, rand);
#endif
#endif //__MULTIRING_FHC_RUNNING__
       
       fxml << "      " << Form("</data_set_%d>\n",iexperiment) << endl;       
     }//save to xml file

     // If a ROOT filename was specified, save the ROOT histogram too
     if(f) {
#ifdef __FHC_RUNNING__
#ifndef __NUEONLY__
       AddExtraEventsToHist(mock_data_1Rmu, full_points_to_add_1Rmu);
#endif
#ifndef __NUMUONLY__
       AddExtraEventsToHist(mock_data_1Re,  full_points_to_add_1Re);
#endif
#endif //__FHC_RUNNING__
       
#ifdef __RHC_RUNNING__	 
#ifndef __NUEONLY__
       AddExtraEventsToHist(mock_data_1Rmu_RHC, full_points_to_add_1Rmu_RHC);
#endif
#ifndef __NUMUONLY__
       AddExtraEventsToHist(mock_data_1Re_RHC,  full_points_to_add_1Re_RHC);
#endif
#endif //__RHC_RUNNING__
       
#ifdef __MULTIRING_FHC_RUNNING__
#ifndef __NUEONLY__
       AddExtraEventsToHist(mock_data_MultiRmu, full_points_to_add_MultiRmu);
#endif
#ifndef __NUMUONLY__
       AddExtraEventsToHist(mock_data_MultiRe,  full_points_to_add_MultiRe);
#endif
#endif //__MULTIRING_FHC_RUNNING__
       
       f->cd();
       
#ifdef __FHC_RUNNING__
#ifndef __NUEONLY__
       mock_data_1Rmu->Write(Form("mock_data_1Rmu_%d",iexperiment));
       delete mock_data_1Rmu;
#endif
#ifndef __NUMUONLY__
       mock_data_1Re->Write(Form("mock_data_1Re_%d",iexperiment));
       delete mock_data_1Re;
#endif
#endif //__FHC_RUNNING__
       
#ifdef __RHC_RUNNING__
#ifndef __NUEONLY__
       mock_data_1Rmu_RHC->Write(Form("mock_data_1Rmu_RHC_%d",iexperiment));
       delete mock_data_1Rmu_RHC;
#endif
#ifndef __NUMUONLY__
       mock_data_1Re_RHC->Write(Form("mock_data_1Re_RHC_%d",iexperiment));
       delete mock_data_1Re_RHC;
#endif
#endif //__RHC_RUNNING__
       
#ifdef __MULTIRING_FHC_RUNNING__
#ifndef __NUEONLY__
       mock_data_MultiRmu->Write(Form("mock_data_MultiRmu_%d",iexperiment));
       delete mock_data_MultiRmu;
#endif
#ifndef __NUMUONLY__
       mock_data_MultiRe->Write(Form("mock_data_MultiRe_%d",iexperiment));
       delete mock_data_MultiRe;
#endif
#endif //__MULTIRING_FHC_RUNNING__

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

  cerr << "Warning h_0, h_1, ... have not been filled (what are they?!)" << endl;
  h_0->Write();
  h_1->Write();
  h_2->Write();
  h_3->Write();
  h_4->Write();
  h_4_Escale->Write();

  //
  // Close files
  //
  if(save_in_xml_format) {
    fxml << "    </observation>       " << endl;
    fxml << "  </superk>              " << endl;
    fxml << "</analysis_inputs>       ";// << endl;
    fxml.close();
  }
  if(f) {
    f->Close();
    delete f;
  }
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
  // Adding Ereco values as a comma separated list.
  // Assuming uniform distribution within each bin, randomize values within each bin.
  for(int ibin = 1; ibin <= h->GetNbinsX(); ibin++) {
    double Nb   = h->GetBinContent (ibin); // number of events in this Ereco bin
    double Emin = h->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
    double Emax = h->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
    double Ecent = h->GetXaxis()->GetBinCenter (ibin); // Ereco bin center
    for(int iev = 0; iev < Nb; iev++) {
      double E= Emin + rand.Uniform() * (Emax-Emin);
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
    fxml << ",\n"<<endl; // should really have an XML comment here but TXMLEngine won't let us.
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
    double Xcent = h->GetXaxis()->GetBinCenter (ibinx); // X bin center
    for(int ibiny = 1; ibiny <= h->GetNbinsY(); ibiny++) {
      double Nb    = h->GetBinContent (ibinx, ibiny); // number of events in this bin
      double Ymin  = h->GetYaxis()->GetBinLowEdge(ibiny); // Y bin min
      double Ymax  = h->GetYaxis()->GetBinUpEdge (ibiny); // Y bin max
      double Ycent = h->GetYaxis()->GetBinCenter (ibiny); // Y bin center
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
    fxml << ",\n" << endl; // should really have an XML comment here but TXMLEngine won't let us.
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
    for(int ip=0; ip < vecp.size(); ip++) {
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

string * ConstructPointsToAddFromXML(TXMLEngine & xml, XMLNodePointer_t & root_node, valor::conv::SuperkSample_t s, string sextra)
{
#ifndef __USE_2D_TEMPLATES__ //1D templates
  string sxml = util::XML2String(&xml, root_node, "superk/observation/data_set_0/one_mu_ring/ereco_values");
  string * sout = new string(sxml);
  if(!sextra.empty()) {
    sout->append(",");
    sout->append(sextra);
  }
  return sout;
#else
  cerr << "TODO - implement ConstructPointsToAddFromXML() in 2D" << endl;
  exit(-1);
#endif
}
