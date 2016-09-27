//
// ------------------------------------------------------------------------------------------------
// Generate / store toy-MC experiments.
// The data can be stored in an XML and/or ROOT format.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

 
// Run in the following way:

//$T2KOSCANA/bin/generate_toymc_ensemble --enable-stat-fluct --syst-tweak-mode 1 --enable-rethrow-unphysical-toys --randomize-all-oscparams --prior-reactors --cpv-nh -3.14 --cpv-ih -3.14 --enable-marg-mass-hierarchy --normal-mass-hierarchy-prior 0.5 --enable-synch-osc-ih-to-nh --file-tag "ciao" --config-file analysis_configuration_test_theta_elike_RFGRPA_Nu6.57E20POT.xml --run 0 --enable-verbose-syst-throws 0 --dmsq-scheme "auto" --mix-param-scheme "sinsqtheta" --cpv-phase-scheme "delta" --mass-hierarchy "normal" --mix-param-octant "two" --hypothesis "3f_osc" --prior 1 --first-experiment 0 --n-experiments 10

#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "T2KPhysicalParameterLimits.h"

#include "Utilities.h"
#include "Binning.h"
#include "AnalysisOptions.h"

#include "Build.h"
#include "UtlForMacrosT2XK.h"

#include "valor/Utilities/UtlForMacros.h"
#include "valor/Oscillation/OscParams.h"
#include "valor/CmdLine/RunOpt_TO_DELETE.h"
#include "valor/Utilities/Env.h"
#include "Conventions/PhysHypothesisId.h"

#include "TMath.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TStopwatch.h"

#include <sstream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace valor;
using namespace valor::conv;

class MyRunOpt : public RunOpt
{
public:
  MyRunOpt(int argc, char ** argv, int parseopt) {
    Init();
    ReadFromCommandLine(argc, argv, parseopt);
    Print(cout);
  }
  ~MyRunOpt() {}


  void   Print (ostream & stream, const int w = 40) const {
    RunOpt::Print(stream, w);
    stream << "Macro-specific running options:" << endl
	   << " - " << setw(w) << "Random seed for position inside bin:" <<  fToySeed << endl
	   << "Output:" << endl
	   << " - " << setw(w) << "Storing toy-MC expts in ROOT format:" <<  ((fSaveAsROOT)     ? "YES" : "NO") << endl
	   << " - " << setw(w) << "Storing toy-MC expts in XML format:"  <<  ((fSaveAsXML)      ? "YES" : "NO") << endl
	   << " - " << setw(w) << "Storing toy-MC fluctuations in XML:"  <<  ((fSaveFluctInXML) ? "YES" : "NO") << endl
	   << "Hybrid data/toy-MC experiments:" << endl
	   << " - " << setw(w) << "Adding extra points from XML:"       << (fXMLExtra.empty() ? "NO" : Form("YES: %s", fXMLExtra.c_str())) << endl
	   << " - " << setw(w) << "Using extra points with dataset ID:" << fXMLExtraId << endl;
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      stream << " - " << setw(w) << Form("Add extra points to %s:", util::AsStringCompact(s).c_str()) << (fExtra[is].empty() ? "NONE" : fExtra[is].c_str()) << endl;
    }
    stream << "Super-K samples to calculate the likelihood weights (penalty term):" << endl;
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      stream << " - " << setw(w) << "Use " << util::AsStringCompact(s) << ": " << (fUseSKsampleWeight[is] ? "YES" : "NO") << endl;
    }
    stream << endl;

  }


  int      ToySeed        () { return fToySeed; }
  bool     SaveAsROOT     () { return fSaveAsROOT; }
  bool     SaveAsXML      () { return fSaveAsXML; }
  bool     SaveFluctInXML () { return fSaveFluctInXML; }
  string   XMLExtraFile   () { return fXMLExtra; }
  int      XMLExtraId     () { return fXMLExtraId; }
  string ExtraPoints(int is) { return fExtra[is]; }
  bool   UseSKsampleWeight(int is){ return fUseSKsampleWeight[is]; }
private:
  int    fToySeed;
  bool   fSaveAsROOT;
  bool   fSaveAsXML;
  bool   fSaveFluctInXML;
  string fXMLExtra;
  int    fXMLExtraId;
  string fExtra[kns];
  bool   fUseSKsampleWeight[kns];
  
  void Init() {
    RunOpt::Init();
    fToySeed        = 0;
    fSaveAsROOT     = true;
    fSaveAsXML      = true;
    fSaveFluctInXML = true;
    fXMLExtra       = "";
    fXMLExtraId     = 0;
    for(int is = 0; is < kns; is++){
      fExtra[is] = "";
      fUseSKsampleWeight[is] = false;
    }
  }
  void ReadFromCommandLine(int argc, char ** argv, int parseopt) {
    RunOpt::ReadFromCommandLine(argc, argv, parseopt);

    if(fParser->OptionExists("toy-random-seed"))
      fToySeed = fParser->ArgAsInt("toy-random-seed");
    if(fParser->OptionExists("enable-save-as-root"))
      fSaveAsROOT = true;
    else if(fParser->OptionExists("disable-save-as-root"))
      fSaveAsROOT = false;
    if(fParser->OptionExists("enable-save-as-xml"))
      fSaveAsXML = true;
    else if(fParser->OptionExists("disable-save-as-xml"))
      fSaveAsXML = false;
    if(fParser->OptionExists("enable-save-fluct-in-xml"))
      fSaveFluctInXML= true;

    else if(fParser->OptionExists("disable-save-fluct-in-xml"))
      fSaveFluctInXML= false;

    if(fParser->OptionExists("xml-file-to-add"))
      fXMLExtra = util::FormFilePath(fParser->ArgAsString("xml-file-to-add").c_str(), "data");
    if(fParser->OptionExists("xml-to-add-dataset-id"))
      fXMLExtraId = fParser->ArgAsInt("xml-to-add-dataset-id");
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      if(fParser->OptionExists(Form("--points-to-add-%s:", util::AsStringCompact(s).c_str())))
	fExtra[is] = fParser->ArgAsString(Form("--points-to-add-%s:", util::AsStringCompact(s).c_str()));
    }//is


    
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      string parse_option = Form("calc-weight-%s", util::AsStringCompact(s).c_str());
      if(fParser->OptionExists(parse_option)){
	fUseSKsampleWeight[is] = true; //
	//cout << util::AsStringCompact(s) << endl;
      }
    }//is

  }
  void Usage(const int w = 40) {
    RunOpt::Usage(w);
    cout << "Macro-specific running options:" << endl
	 << " " << setw(w) << "--toy-random-seed"           << "Random seed used to specify where inside the bin the toys are generated" << endl
	 << "Output:" << endl
	 << " " << setw(w) << "--enable-save-as-root"       << "Save in ROOT format?" << endl
	 << " " << setw(w) << "--enable-save-as-xml"        << "Save in XML format?" << endl
	 << " " << setw(w) << "--enable-save-fluct-in-xml"  << "Save fluctuations in the XML file?" << endl
	 << "Hybrid data/toy-MC experiments:" << endl
	 << " " << setw(w) << "--xml-file-to-add"        << "XML file (in $T2KOSCANA/inputs/) with data to add" << endl
	 << " " << setw(w) << "--xml-to-add-dataset-id"  << "Dataset ID in the XML file to add" << endl;
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      cout << " " << setw(w) << Form("--points-to-add-%s:", util::AsStringCompact(s).c_str()) << "List of data-points to add eg. E1:q1,E2:q2,E3:q3" << endl;
    }

    cout << "Super-K samples to calculate the likelihood weights (penalty term):" << endl;
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      cout << " - " << setw(w) << "Use " << util::AsStringCompact(s) << ": " << (fUseSKsampleWeight[is] ? "YES" : "NO") << endl;
    }

    cout << endl;
  }
};


double calc_chisq(TH2D* hobs, TH2D* hexp){
  int nbins_x = hobs->GetNbinsX();
  int nbins_y = hobs->GetNbinsY();
  double chisq = 0;
  
  for(int i=1;i<= nbins_x;i++)
    {
      for(int j=1;j<= nbins_y;j++)
	{
	  double obs = hobs->GetBinContent(i,j);
	  double exp = hexp->GetBinContent(i,j);
	  chisq += valor::util::ChisqPoisson(obs,exp);
	}
    } 
  return chisq;
}


int main(int argc, char ** argv)
{
  TStopwatch sw;

  Env::Instance()->AutoSet("t2xk");
  cout << *Env::Instance() << endl;
  gROOT->SetBatch(true);
  MyRunOpt runopt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParse3FOscParams );

  //
  
  TRandom3 rand(runopt.ToySeed()); //random number to see where inside a bin a toy is generated

  // Extra points to add to the toy dataset(s)
  string * full_points_to_add[kns] = { 0,0,0,0,0,0 };
  // Are we adding events from an XML file?
  if (!runopt.XMLExtraFile().empty()) {
    TXMLEngine xml;
    XMLDocPointer_t xmldoc = xml.ParseFile(runopt.XMLExtraFile().c_str());
    if (xmldoc == 0) {
      cerr << "Can not parse XML file to add: " << runopt.XMLExtraFile() << endl;
      exit(1);
    }
    XMLNodePointer_t root_node = xml.DocGetRootElement(xmldoc);
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      full_points_to_add[is] = valor::util::ConstructPointsToAddFromXML(xml, root_node, s, runopt.ExtraPoints(is), runopt.XMLExtraId());
    }//is
  }//runopt.XMLExtraFile()
  else {
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      if ( ! runopt.ExtraPoints(is).empty() ) full_points_to_add[is] = new string(runopt.ExtraPoints(is));
    }//is
  }//!runopt.XMLExtraFile()
  cout << "   - Adding the following data-points in each toy-MC expt: " << endl;
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    if (s == conv::kSKsUndefined) continue;
    cout << "    For " << util::AsStringCompact(s) << "\t" << ((full_points_to_add[is]) ? *full_points_to_add[is] : "*None*") << endl;
  }//is
  cout << endl;
  
  // Open files
  TFile * f = 0;
  if(runopt.SaveAsROOT()) {
    ostringstream root_filename;
    root_filename << runopt.FileTag() << "-" << runopt.Run() << ".root";
    cout << root_filename.str();
    f = new TFile(root_filename.str().c_str(), "recreate");
  }

  ofstream fxml;
  ofstream ftxt;
  if(runopt.SaveAsXML()) {
    ostringstream xml_filename;
    xml_filename << runopt.FileTag() << "-" << runopt.Run() << ".xml";
    fxml.open(xml_filename.str().c_str());
    fxml << "<?xml version=\"1.0\"?> " << endl;
    fxml << "<analysis_inputs>       " << endl;
    fxml << "  <superk>              " << endl;
    fxml << "    <observation>       " << endl;
    ostringstream txt_filename;
    txt_filename << runopt.FileTag() << "-" << runopt.Run() << ".txt";
    ftxt.open(txt_filename.str().c_str());
  }

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(runopt.ConfigFile().c_str());



  // Initialize likelihood weight histogram (used as penalty term) 
  
  TH1D *hChisqWeight = new TH1D("hChisqWeight","hChisqWeight",4000,0,400);
  //TH1D *hLikeWeight = new TH1D("hLikeWeight","hLikeWeight",4000,0,0.000000000001);
  
  // Data files
  
  TH2D *h_data_1rmu;
  TH2D *h_data_1re;
  TH2D *h_data_rhc1rmu;
  TH2D *h_data_rhc1re;
  
  string rootfile = runopt.PathFileData();
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile.c_str()));
  SuperkObservation * skobs_fakedata = new SuperkObservation;
  TFile *infile;
  cout << endl;
  cout << " - use root file" << rootfile2 << endl;
  cout << endl;
  infile = new TFile(rootfile2.Data(),"READ");
  if(!infile->IsOpen()){
    cout << endl;
    cout << "The data root file doesn't exist: don't draw data spectra" << endl;
    cout << endl;
    infile->Close();
  }
  conv::SuperkSample_t s_1rmu    = util::SKSampleFromString("1Rmu"); 
  conv::SuperkSample_t s_1re     = util::SKSampleFromString("1Re"); 
  conv::SuperkSample_t s_rhc1rmu = util::SKSampleFromString("RHC1Rmu"); 
  conv::SuperkSample_t s_rhc1re  = util::SKSampleFromString("RHC1Re"); 
  
  //int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  bool using_xml_file = false;
  int exp_id = 0;
  skobs_fakedata->Read(rootfile2.Data(),exp_id,using_xml_file);
  h_data_1rmu    = (TH2D*) skobs_fakedata->TwoDSpectrum(s_1rmu)   ->Clone("h_data_1rmu");
  h_data_1re     = (TH2D*) skobs_fakedata->TwoDSpectrum(s_1re)    ->Clone("h_data_1re");
  h_data_rhc1rmu = (TH2D*) skobs_fakedata->TwoDSpectrum(s_rhc1rmu)->Clone("h_data_rhc1rmu");
  h_data_rhc1re  = (TH2D*) skobs_fakedata->TwoDSpectrum(s_rhc1re) ->Clone("h_data_rhc1re"); 
  
  cout << endl;
  cout << "Take data set from: " << runopt.PathFileData() << endl;
  cout << "# of events of Super-K samples:" << endl;
  cout << " - 1Rmu: "    << h_data_1rmu->Integral() << endl;
  cout << " - 1Re: "     << h_data_1re->Integral() << endl;
  cout << " - RHC1Rmu: " << h_data_rhc1rmu->Integral() << endl;
  cout << " - RHC1Re: "  << h_data_rhc1re->Integral() << endl;
  cout << endl;
  
  //
 

  int weird_applied_neg      = 0;
  int weird_unphysical_tweak = 0;
  
  //for(int iexperiment = runopt.FirstExperiment(), last_experiment = runopt.LastExperiment(); iexperiment <= last_experiment; iexperiment++) {
  
  int iexperiment = runopt.FirstExperiment();
  int last_experiment = runopt.LastExperiment();
  while( iexperiment <= last_experiment ){
    
    cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;

     // Generate mock data
     cout << "- Generating toy-MC..." << endl;
     SuperkObservation * skobs = new SuperkObservation;
     PhysHypothesis ph(runopt.Hypothesis());


     //util::OscParams::Instance()->Print(cout);

     // Randomize specified osc params for current toy MC expt
     util::OscParams::Instance()->Randomize();

     // Set the mass hierarchy
     util::SetMassHierarchy(skexp);

     // Update physics hypothesis with chosen value of oscillation parameters
     util::OscParams::Instance()->UpdatePhysHypothesis(ph);

     // Fill toy histograms of osc params
     util::OscParams::Instance()->FillOscToyHisto(ph);

     bool applied_negevents_condition = true;
     bool unphysical_tweak = true;
     while(applied_negevents_condition || unphysical_tweak) {

       util::PerformSysThrow(runopt.SystTweakMode(), skexp, runopt.VerboseSystThrows());

       // Rethrow unphysical parameters
       PhysParamList true_nuisance_params = skexp->Tweaks();
       if (runopt.RethrowIfUnphysical()) {
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

       applied_negevents_condition = skobs->CreateMockData(*skexp, ph, runopt.StatFluctuations());
       cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
       if(applied_negevents_condition) weird_applied_neg++;
     }//applied_negevents_condition


     //Get the spectra
#ifndef __USE_2D_TEMPLATES__ //1D templates
     TH1D * mock_data[kns];
     TH1D * mock_exp_data[kns];
     TH1D * data_real[kns];
#else
     TH2D * mock_data[kns];
     TH2D * mock_exp [kns];
     TH2D * data_real[kns];
#endif
    //Get exp data
#ifndef __USE_2D_TEMPLATES__ //1D templates
    TH1D * exp_data[kns];
#else
    TH2D * exp_data[kns];
#endif
     for (int is = 0 ; is < conv::kns ; is++) {
       SuperkSample_t s = conv::kSKSamples[is];
       if (s == conv::kSKsUndefined) continue;
#ifndef __USE_2D_TEMPLATES__ //1D templates
       mock_data[is] = skobs->ErecoSpectrum(s);
       mock_exp[is]  = skexp->CurrPredictedErecoSpectrum(s);
#else
       mock_data[is] = skobs->TwoDSpectrum(s);
       mock_exp[is]  = skexp->CurrPredictedErecoVsQ2Spectrum(s);
#endif
#ifndef __USE_2D_TEMPLATES__
       exp_data[is] = skexp->CurrPredictedErecoSpectrum(s);
#else
       exp_data[is] = skexp->CurrPredictedErecoVsQ2Spectrum(s);
#endif
       exp_data[is]->SetDirectory(0);
       mock_data[is]->SetDirectory(0);
       mock_exp[is] ->SetDirectory(0);
     }
     
     // Get the data
     data_real[s_1rmu]    = h_data_1rmu;
     data_real[s_1re]     = h_data_1re;
     data_real[s_rhc1rmu] = h_data_rhc1rmu;
     data_real[s_rhc1re]  = h_data_rhc1re;

     
     // Calculate likelihood weights based on the real data set
  
     double toy_weight = 1.;
     double toy_chisq = 0.;
     
     for (int is = 0 ; is < conv::kns ; is++) {
       SuperkSample_t s = conv::kSKSamples[is];
       if (s == conv::kSKsUndefined) continue;


       if(runopt.UseSKsampleWeight(is)){
	 toy_chisq += calc_chisq(data_real[is],mock_exp[is]);    
	 
	 // cout << util::AsString(s) << ": " 
	 //      << data_real[is]->Integral() << " - "
	 //      << mock_data[is]->Integral() << " - "
	 //      << mock_exp[is]->Integral() << endl;
	 // cout << "chisq = " << calc_chisq(data_real[is],mock_exp[is]) << endl;
       }
     }
     toy_weight = TMath::Exp(-1*toy_chisq/2);

     double unif_weight = rand.Uniform(0,1);
     
     if(toy_weight<unif_weight){
       
       cout << "like_unif = " << unif_weight << " > " << "like_weight = " << toy_weight << " --> reject toy" << endl;
 
       continue; // if like_fakedata < like_rndm --> reject the toy
     }

     //cout << "toy_chisq = " << toy_chisq << endl;
     //cout << "toy_weight = " << toy_weight << endl;

     hChisqWeight->Fill(toy_chisq);
     
     // If an XML file was specified, add an entry
     if(runopt.SaveAsXML()) {

       // Dataset info
       fxml << "      " << Form("<data_set_%d>",iexperiment) << endl;
       if(runopt.SaveFluctInXML()) {
         fxml << "<!-- ";
         ph.Print(fxml);
         if(runopt.SystTweakMode())
           skexp->PrintSystematicVariations(fxml);
         fxml << " -->" << endl;
       }

       // The 'observed' data
       for (int is = 0 ; is < conv::kns ; is++) {
	 SuperkSample_t s = conv::kSKSamples[is];
	 if (s == conv::kSKsUndefined) continue;
	 valor::util::AddEventsToXML(fxml, ftxt, mock_data[is], s, full_points_to_add[is], rand);
       }//is
       
       fxml << "      " << Form("</data_set_%d>\n",iexperiment) << endl;       
     }//save to xml file

     // If a ROOT filename was specified, save the ROOT histogram too
     if(f) {
       for (int is = 0 ; is < conv::kns ; is++) {
	 SuperkSample_t s = conv::kSKSamples[is];
	 if (s == conv::kSKsUndefined) continue;
	 valor::util::AddExtraEventsToHist(mock_data[is], full_points_to_add[is]);
       }//is
      
       f->cd();

       for (int is = 0 ; is < conv::kns ; is++) {
         SuperkSample_t s = conv::kSKSamples[is];
         if (s == conv::kSKsUndefined) continue;
         mock_data[is]->Write(Form("mock_data_%s_%d", util::AsStringCompact(s).c_str(), iexperiment));
         delete mock_data[is];
         exp_data[is]->Write(Form("exp_data_%s_%d", util::AsStringCompact(s).c_str(), iexperiment));
         delete exp_data[is];
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

     iexperiment++;

  }//iexperiment


  hChisqWeight->Write();


  //
  // Close files
  //
  if(runopt.SaveAsXML()) {
    fxml << "    </observation>       " << endl;
    fxml << "  </superk>              " << endl;
    fxml << "</analysis_inputs>       ";// << endl;
    fxml.close();
  }
  if(f) {
    util::OscParams::Instance()->Write(f);
    f->Close();
    delete f;
  }//f

  cout << endl 
       << " - Summary - " << endl
       << "    * Created " << runopt.NExperiments() << " toys" << endl
       << "    * " << weird_unphysical_tweak << " unphysical_tweaks rethrown" << endl
       << "    * " << weird_applied_neg << " applied_negevents_condition rethrown" << endl
       << endl;

  sw.Stop();
  sw.Print();
}



