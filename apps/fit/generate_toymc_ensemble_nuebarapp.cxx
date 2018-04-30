//
// ------------------------------------------------------------------------------------------------
// Generate / store toy-MC experiments.
// The data can be stored in an XML and/or ROOT format.
// This is a special script which can generate toys based on priors but also takes a data rootfile
// the data rootfile is used for the 3 samples 1Rmu, 1Re, RHC1Rmu. The expectations and the data
// for the randomly thrown toy are saved along with the data from the 3 samples.

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

#include "valor/DerivedUtilities/UtlForMacros.h"
#include "valor/Oscillation/OscParams.h"
#include "valor/CmdLine/RunOpt_TO_DELETE.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Conventions/PhysHypothesisId.h"

#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TStopwatch.h"

#include <sstream>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;
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
    stream << endl;
  }
  int      ToySeed        () { return fToySeed; }
  bool     SaveAsROOT     () { return fSaveAsROOT; }
  bool     SaveAsXML      () { return fSaveAsXML; }
  bool     SaveFluctInXML () { return fSaveFluctInXML; }
  string   XMLExtraFile   () { return fXMLExtra; }
  int      XMLExtraId     () { return fXMLExtraId; }
  string ExtraPoints(int is) { return fExtra[is]; }
private:
  int    fToySeed;
  bool   fSaveAsROOT;
  bool   fSaveAsXML;
  bool   fSaveFluctInXML;
  string fXMLExtra;
  int    fXMLExtraId;
  string fExtra[kns];
  
  void Init() {
    RunOpt::Init();
    fToySeed        = 0;
    fSaveAsROOT     = true;
    fSaveAsXML      = true;
    fSaveFluctInXML = true;
    fXMLExtra       = "";
    fXMLExtraId     = 0;
    for(int is = 0; is < kns; is++)
      fExtra[is] = "";
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
    cout << endl;
  }
};

//____________________________________________________________________________
TH1D *hToyHisto_dcp_nh;
TH1D *hToyHisto_Mix13_nh;
TH1D *hToyHisto_Mix23_nh;
TH1D *hToyHisto_Dm23_nh;
TH1D *hToyHisto_Mix12_nh;
TH1D *hToyHisto_Dm12_nh;
TH1D *hToyHisto_dcp_ih;
TH1D *hToyHisto_Mix13_ih;
TH1D *hToyHisto_Mix23_ih;
TH1D *hToyHisto_Dm23_ih;
TH1D *hToyHisto_Mix12_ih;
TH1D *hToyHisto_Dm12_ih;
TH1D *hToyHisto_MH;

//NuBar
TH1D *hToyHisto_dcpBar_nh;
TH1D *hToyHisto_Mix13Bar_nh;
TH1D *hToyHisto_Mix23Bar_nh;
TH1D *hToyHisto_Dm23Bar_nh;
TH1D *hToyHisto_Mix12Bar_nh;
TH1D *hToyHisto_Dm12Bar_nh;
TH1D *hToyHisto_dcpBar_ih;
TH1D *hToyHisto_Mix13Bar_ih;
TH1D *hToyHisto_Mix23Bar_ih;
TH1D *hToyHisto_Dm23Bar_ih;
TH1D *hToyHisto_Mix12Bar_ih;
TH1D *hToyHisto_Dm12Bar_ih;

//____________________________________________________________________________
// Re-implemented from valor::util::OscParams::FillOscToyHisto (as that class is a singleton and I need to look at rejected throws too)
// TODO: Do this in a smarter way instead of re-implementing the needed methods
void FillOscToyHisto(PhysHypothesis &ph)
{
  // Fill toy osc params histograms
  double mix12 = conv::kBadResult;
  double mix13 = conv::kBadResult;
  double mix23 = conv::kBadResult;
  double dm12  = conv::kBadResult;
  double dm23  = conv::kBadResult;
  double dcp   = conv::kBadResult;
  double mix12bar = conv::kBadResult;
  double mix13bar = conv::kBadResult;
  double mix23bar = conv::kBadResult;
  double dm12bar  = conv::kBadResult;
  double dm23bar  = conv::kBadResult;
  double dcpbar   = conv::kBadResult;

  if(ph.Contains(kPhypMixing12)) mix12 = ph.CurrValue(kPhypMixing12);
  if(ph.Contains(kPhypMixing13)) mix13 = ph.CurrValue(kPhypMixing13);
  if(ph.Contains(kPhypMixing23)) mix23 = ph.CurrValue(kPhypMixing23);
  if(ph.Contains(kPhypDm12sq))   dm12 = ph.CurrValue(kPhypDm12sq);
  if(ph.Contains(kPhypDm23sq))   dm23 = ph.CurrValue(kPhypDm23sq);
  if(ph.Contains(kPhypCPV))      dcp = ph.CurrValue(kPhypCPV);

  if(ph.Contains(kPhypMixing12Bar)) mix12bar = ph.CurrValue(kPhypMixing12Bar);
  if(ph.Contains(kPhypMixing13Bar)) mix13bar = ph.CurrValue(kPhypMixing13Bar);
  if(ph.Contains(kPhypMixing23Bar)) mix23bar = ph.CurrValue(kPhypMixing23Bar);
  if(ph.Contains(kPhypDm12sqBar))   dm12bar = ph.CurrValue(kPhypDm12sqBar);
  if(ph.Contains(kPhypDm23sqBar))   dm23bar = ph.CurrValue(kPhypDm23sqBar);
  if(ph.Contains(kPhypCPVBar))      dcpbar = ph.CurrValue(kPhypCPVBar);

  int mh = Opt::Instance()->gMassHierarchyAssumed;

  if(mh==conv::kNmhNormal){
    hToyHisto_Mix12_nh->Fill(mix12);
    hToyHisto_Mix13_nh->Fill(mix13);
    hToyHisto_Mix23_nh->Fill(mix23);
    hToyHisto_Dm12_nh->Fill(dm12);
    hToyHisto_Dm23_nh->Fill(dm23);
    hToyHisto_dcp_nh->Fill(dcp);
    hToyHisto_Mix12Bar_nh->Fill(mix12bar);
    hToyHisto_Mix13Bar_nh->Fill(mix13bar);
    hToyHisto_Mix23Bar_nh->Fill(mix23bar);
    hToyHisto_Dm12Bar_nh->Fill(dm12bar);
    hToyHisto_Dm23Bar_nh->Fill(dm23bar);
    hToyHisto_dcpBar_nh->Fill(dcpbar);
  }
  else if(mh==conv::kNmhInverted){
    hToyHisto_Mix12_ih->Fill(mix12);
    hToyHisto_Mix13_ih->Fill(mix13);
    hToyHisto_Mix23_ih->Fill(mix23);
    hToyHisto_Dm12_ih->Fill(dm12);
    hToyHisto_Dm23_ih->Fill(dm23);
    hToyHisto_dcp_ih->Fill(dcp);
    hToyHisto_Mix12Bar_ih->Fill(mix12bar);
    hToyHisto_Mix13Bar_ih->Fill(mix13bar);
    hToyHisto_Mix23Bar_ih->Fill(mix23bar);
    hToyHisto_Dm12Bar_ih->Fill(dm12bar);
    hToyHisto_Dm23Bar_ih->Fill(dm23bar);
    hToyHisto_dcpBar_ih->Fill(dcpbar);
  }
  hToyHisto_MH->Fill(mh);
}
//____________________________________________________________________________
// Re-implemented from valor::util::OscParams::Write
void WriteAllThrows(TFile * f)
{
  f->cd();

  //toys used
  if(hToyHisto_dcp_nh) {
    hToyHisto_dcp_nh->SetDirectory(0);
    hToyHisto_dcp_nh->Write();
  }
  if(hToyHisto_Mix13_nh) {
    hToyHisto_Mix13_nh->SetDirectory(0);
    hToyHisto_Mix13_nh->Write();
  }
  if(hToyHisto_Mix23_nh) {
    hToyHisto_Mix23_nh->SetDirectory(0);
    hToyHisto_Mix23_nh->Write();
  }
  if(hToyHisto_Dm23_nh) {
    hToyHisto_Dm23_nh->SetDirectory(0);
    hToyHisto_Dm23_nh->Write();
  }
  if(hToyHisto_Mix12_nh) {
    hToyHisto_Mix12_nh->SetDirectory(0);
    hToyHisto_Mix12_nh->Write();
  }
  if(hToyHisto_Dm12_nh) {
    hToyHisto_Dm12_nh->SetDirectory(0);
    hToyHisto_Dm12_nh->Write();
  }
  if(hToyHisto_dcp_ih) {
    hToyHisto_dcp_ih->SetDirectory(0);
    hToyHisto_dcp_ih->Write();
  }
  if(hToyHisto_Mix13_ih) {
    hToyHisto_Mix13_ih->SetDirectory(0);
    hToyHisto_Mix13_ih->Write();
  }
  if(hToyHisto_Mix23_ih) {
    hToyHisto_Mix23_ih->SetDirectory(0);
    hToyHisto_Mix23_ih->Write();
  }
  if(hToyHisto_Dm23_ih) {
    hToyHisto_Dm23_ih->SetDirectory(0);
    hToyHisto_Dm23_ih->Write();
  }
  if(hToyHisto_Mix12_ih) {
    hToyHisto_Mix12_ih->SetDirectory(0);
    hToyHisto_Mix12_ih->Write();
  }
  if(hToyHisto_Dm12_ih) {
    hToyHisto_Dm12_ih->SetDirectory(0);
    hToyHisto_Dm12_ih->Write();
  }
  if(hToyHisto_MH) {
    hToyHisto_MH->SetDirectory(0);
    hToyHisto_MH->Write();
  }
  if(hToyHisto_dcpBar_nh) {
    hToyHisto_dcpBar_nh->SetDirectory(0);
    hToyHisto_dcpBar_nh->Write();
  }
  if(hToyHisto_Mix13Bar_nh) {
    hToyHisto_Mix13Bar_nh->SetDirectory(0);
    hToyHisto_Mix13Bar_nh->Write();
  }
  if(hToyHisto_Mix23Bar_nh) {
    hToyHisto_Mix23Bar_nh->SetDirectory(0);
    hToyHisto_Mix23Bar_nh->Write();
  }
  if(hToyHisto_Dm23Bar_nh) {
    hToyHisto_Dm23Bar_nh->SetDirectory(0);
    hToyHisto_Dm23Bar_nh->Write();
  }
  if(hToyHisto_Mix12Bar_nh) {
    hToyHisto_Mix12Bar_nh->SetDirectory(0);
    hToyHisto_Mix12Bar_nh->Write();
  }
  if(hToyHisto_Dm12Bar_nh) {
    hToyHisto_Dm12Bar_nh->SetDirectory(0);
    hToyHisto_Dm12Bar_nh->Write();
  }
  if(hToyHisto_dcpBar_ih) {
    hToyHisto_dcpBar_ih->SetDirectory(0);
    hToyHisto_dcpBar_ih->Write();
  }
  if(hToyHisto_Mix13Bar_ih) {
    hToyHisto_Mix13Bar_ih->SetDirectory(0);
    hToyHisto_Mix13Bar_ih->Write();
  }
  if(hToyHisto_Mix23Bar_ih) {
    hToyHisto_Mix23Bar_ih->SetDirectory(0);
    hToyHisto_Mix23Bar_ih->Write();
  }
  if(hToyHisto_Dm23Bar_ih) {
    hToyHisto_Dm23Bar_ih->SetDirectory(0);
    hToyHisto_Dm23Bar_ih->Write();
  }
  if(hToyHisto_Mix12Bar_ih) {
    hToyHisto_Mix12Bar_ih->SetDirectory(0);
    hToyHisto_Mix12Bar_ih->Write();
  }
  if(hToyHisto_Dm12Bar_ih) {
    hToyHisto_Dm12Bar_ih->SetDirectory(0);
    hToyHisto_Dm12Bar_ih->Write();
  }

}
//____________________________________________________________________________



int main(int argc, char ** argv)
{
  TStopwatch sw;

  Env::Instance()->AutoSet("t2xk");
  cout << *Env::Instance() << endl;
  gROOT->SetBatch(true);
  MyRunOpt *runopt = new MyRunOpt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParseFitMarg | kParse3FOscParams);
  if(runopt->CountFailedToys()){
    hToyHisto_dcp_nh = new TH1D("hToyHisto_dcp_nh_accepted","hToyHisto_dcp_nh_accepted",100,-3.1416,+3.1416);
    hToyHisto_Mix13_nh = new TH1D("hToyHisto_Mix13_nh_accepted","hToyHisto_Mix13_nh_accepted",100,0.,0.04);
    hToyHisto_Mix23_nh = new TH1D("hToyHisto_Mix23_nh_accepted","hToyHisto_Mix23_nh_accepted",100,0.3,0.7);
    hToyHisto_Dm23_nh = new TH1D("hToyHisto_Dm23_nh_accepted","hToyHisto_Dm23_nh_accepted",100,2.3E-3,2.7E-3);
    hToyHisto_Mix12_nh = new TH1D("hToyHisto_Mix12_nh_accepted","hToyHisto_Mix12_nh_accepted",100,0.2,0.4);
    hToyHisto_Dm12_nh = new TH1D("hToyHisto_Dm12_nh_accepted","hToyHisto_Dm12_nh_accepted",100,6.E-5,9.E-5);
    hToyHisto_dcp_ih = new TH1D("hToyHisto_dcp_ih_accepted","hToyHisto_dcp_ih_accepted",100,-3.1416,+3.1416);
    hToyHisto_Mix13_ih = new TH1D("hToyHisto_Mix13_ih_accepted","hToyHisto_Mix13_ih_accepted",100,0.,0.04);
    hToyHisto_Mix23_ih = new TH1D("hToyHisto_Mix23_ih_accepted","hToyHisto_Mix23_ih_accepted",100,0.3,0.7);
    hToyHisto_Dm23_ih = new TH1D("hToyHisto_Dm23_ih_accepted","hToyHisto_Dm23_ih_accepted",100,2.3E-3,2.7E-3);
    hToyHisto_Mix12_ih = new TH1D("hToyHisto_Mix12_ih_accepted","hToyHisto_Mix12_ih_accepted",100,0.2,0.4);
    hToyHisto_Dm12_ih = new TH1D("hToyHisto_Dm12_ih_accepted","hToyHisto_Dm12_ih_accepted",100,6.E-5,9.E-5);
    hToyHisto_dcpBar_nh = new TH1D("hToyHisto_dcpBar_nh_accepted","hToyHisto_dcpBar_nh_accepted",100,-3.1416,+3.1416);
    hToyHisto_Mix13Bar_nh = new TH1D("hToyHisto_Mix13Bar_nh_accepted","hToyHisto_Mix13Bar_nh_accepted",100,0.,0.04);
    hToyHisto_Mix23Bar_nh = new TH1D("hToyHisto_Mix23Bar_nh_accepted","hToyHisto_Mix23Bar_nh_accepted",100,0.3,0.7);
    hToyHisto_Dm23Bar_nh = new TH1D("hToyHisto_Dm23Bar_nh_accepted","hToyHisto_Dm23Bar_nh_accepted",100,2.3E-3,2.7E-3);
    hToyHisto_Mix12Bar_nh = new TH1D("hToyHisto_Mix12Bar_nh_accepted","hToyHisto_Mix12Bar_nh_accepted",100,0.2,0.4);
    hToyHisto_Dm12Bar_nh = new TH1D("hToyHisto_Dm12Bar_nh_accepted","hToyHisto_Dm12Bar_nh_accepted",100,6.E-5,9.E-5);
    hToyHisto_dcpBar_ih = new TH1D("hToyHisto_dcpBar_ih_accepted","hToyHisto_dcpBar_ih_accepted",100,-3.1416,+3.1416);
    hToyHisto_Mix13Bar_ih = new TH1D("hToyHisto_Mix13Bar_ih_accepted","hToyHisto_Mix13Bar_ih_accepted",100,0.,0.04);
    hToyHisto_Mix23Bar_ih = new TH1D("hToyHisto_Mix23Bar_ih_accepted","hToyHisto_Mix23Bar_ih_accepted",100,0.3,0.7);
    hToyHisto_Dm23Bar_ih = new TH1D("hToyHisto_Dm23Bar_ih_accepted","hToyHisto_Dm23Bar_ih_accepted",100,2.3E-3,2.7E-3);
    hToyHisto_Mix12Bar_ih = new TH1D("hToyHisto_Mix12Bar_ih_accepted","hToyHisto_Mix12Bar_ih_accepted",100,0.2,0.4);
    hToyHisto_Dm12Bar_ih = new TH1D("hToyHisto_Dm12Bar_ih_accepted","hToyHisto_Dm12Bar_ih_accepted",100,6.E-5,9.E-5);
    hToyHisto_MH = new TH1D("hToyHisto_MH_accepted","hToyHisto_MH_accepted",2,0,2);
  }
  
  double scale_toys_based_on_samples_scale =  runopt->Scale_toys_based_on_samples_scale();
  TH1D* chisq_distribution          = new TH1D("chisq_distribution"         ,"chisq_distribution"         ,100000,0,1000);
  TH1D* chisq_distribution_accepted = new TH1D("chisq_distribution_accepted","chisq_distribution_accepted",100000,0,1000);

  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), runopt->ConfigFile().c_str() ));
  TString xmlfile2(runopt->Xmlfile().c_str() );
  TString rootfile2(runopt->Rootfile().c_str());
  
  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(runopt->GlobalMassHierarchy());
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(runopt->GlobalDmsqScheme());
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(runopt->GlobalMixParamScheme());
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(runopt->GlobalMixParamOctant());
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  Opt::Instance()->SetDefaultParameterValues();
  
  if(!runopt->DoToys()){
    cout << "- Reading analysis configuration from: " << config << endl;
    cout << "- Toy-MC" << endl;
    if(runopt->Using_xml_file()) cout << "   - Reading experiments from: " << xmlfile2 << endl;
    else cout << "   - Reading experiments from: " << rootfile2 << endl;
    cout << "   - Using experiments in range: [" << 0 << ", " << 0 << "]" << endl;
  }
  
  // Extra points to add to the toy dataset(s)
  string * full_points_to_add[kns] = { 0,0,0,0,0,0 };
  // Are we adding events from an XML file?
  if (!runopt->XMLExtraFile().empty()) {
    TXMLEngine xml;
    XMLDocPointer_t xmldoc = xml.ParseFile(runopt->XMLExtraFile().c_str());
    if (xmldoc == 0) {
      cerr << "Can not parse XML file to add: " << runopt->XMLExtraFile() << endl;
      exit(1);
    }
    XMLNodePointer_t root_node = xml.DocGetRootElement(xmldoc);
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      full_points_to_add[is] = valor::util::ConstructPointsToAddFromXML(xml, root_node, s, runopt->ExtraPoints(is), runopt->XMLExtraId());
    }//is
  }//runopt->XMLExtraFile()
  else {
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      if ( ! runopt->ExtraPoints(is).empty() ) full_points_to_add[is] = new string(runopt->ExtraPoints(is));
    }//is
  }//!runopt->XMLExtraFile()
  cout << "   - Adding the following data-points in each toy-MC expt: " << endl;
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    if (s == conv::kSKsUndefined) continue;
    cout << "    For " << util::AsStringCompact(s) << "\t" << ((full_points_to_add[is]) ? *full_points_to_add[is] : "*None*") << endl;
  }//is
  cout << endl;
  
  // Open files
  TFile * f = 0;
  if(runopt->SaveAsROOT()) {
    ostringstream root_filename;
    root_filename << runopt->FileTag() << "-" << runopt->Run() << ".root";
    cout << root_filename.str();
    f = new TFile(root_filename.str().c_str(), "recreate");
  }
  
  ofstream fxml;
  ofstream ftxt;
  if(runopt->SaveAsXML()) {
    ostringstream xml_filename;
    xml_filename << runopt->FileTag() << "-" << runopt->Run() << ".xml";
    fxml.open(xml_filename.str().c_str());
    fxml << "<?xml version=\"1.0\"?> " << endl;
    fxml << "<analysis_inputs>       " << endl;
    fxml << "  <superk>              " << endl;
    fxml << "    <observation>       " << endl;
    ostringstream txt_filename;
    txt_filename << runopt->FileTag() << "-" << runopt->Run() << ".txt";
    ftxt.open(txt_filename.str().c_str());
  }
  
  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(runopt->ConfigFile().c_str());
  
  SuperkObservation * skobs_data = new SuperkObservation();
  
  if(!runopt->DoToys()){
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined || s == conv::kSKsRHCFCFV1RingELike) continue;
      
      if(runopt->Using_xml_file()){
        cout << "Use xml file" << endl; // Only read 1st toy
        skobs_data->ReadSample(xmlfile2.Data(),0,runopt->Using_xml_file(),s);
      }
      else{
        cout << "Use root file" << endl; // Only read 1st toy
        skobs_data->ReadSample(rootfile2.Data(),0,runopt->Using_xml_file(),s);
      }
    }
  }
  
  //Compare input data file with each expectation toy to weight the fake data
  TH2D * real_data[kns];
  
  if(runopt->Scale_toys_based_on_samples()){
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined || s == conv::kSKsRHCFCFV1RingELike) continue;
      real_data[is] = skobs_data->TwoDSpectrum(s);
    }
  }
  
  int scale_throws_accepted   = 0;
  int scale_throws_rejected   = 0;
  int weird_applied_neg       = 0;
  int weird_unphysical_tweak  = 0;
  
  for(int iexperiment = runopt->FirstExperiment(), last_experiment = runopt->LastExperiment(); iexperiment <= last_experiment; iexperiment++) {
    cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;
    
    // Set Seeds
    // Because the seed is dependant on iexperiment, using the CountFailedToys option will cause the generated toys to
    // be different than if this option was not used.
    Rnd* rand = Rnd::Instance();
    if(runopt->ToySeed() != 0){
      // 0 is a special case - set all Rnd::Instance() according to the input seed and iexperiment unless the input seed is 0
      UInt_t seed = runopt->ToySeed() + iexperiment;
      rand->SetSeeds(seed, seed+1, seed+2, seed+3, seed+4);
    }  

    // Generate mock data
    cout << "- Generating toy-MC..." << endl;
    SuperkObservation * skobs = new SuperkObservation();
    PhysHypothesis ph(runopt->Hypothesis());
    
    bool applied_negevents_condition = true;
    bool unphysical_tweak = true;
    while(applied_negevents_condition || unphysical_tweak) {
      
      //Randomise oscillaiton parameters must be done inside this loop so
      //they are chosen wrt to MC sampling method.
      // Randomize specified osc params for current toy MC expt
      util::OscParams::Instance()->Randomize();
      
      // Set the mass hierarchy
      util::SetMassHierarchy(skexp);
      
      // Update physics hypothesis with chosen value of oscillation parameters
      util::OscParams::Instance()->UpdatePhysHypothesis(ph);
      
      util::PerformSysThrow(runopt->SystTweakMode(), skexp, runopt->VerboseSystThrows());
      
      // Rethrow unphysical parameters
      PhysParamList true_nuisance_params = skexp->Tweaks();
      if (runopt->RethrowIfUnphysical()) {
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
      
      applied_negevents_condition = skobs->CreateMockData(*skexp, ph, runopt->StatFluctuations());
      cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
      if(applied_negevents_condition)
        weird_applied_neg++;
      
      bool use_it = 1;
      
      //Calc chisq of data and expectation
      //Take toy based on chisq and random number
      if(runopt->Scale_toys_based_on_samples()){
        //Get exp data
#ifndef __USE_2D_TEMPLATES__ //1D templates
        TH1D * exp_data[kns];
#else
        TH2D * exp_data[kns];
#endif
        for (int is = 0 ; is < conv::kns ; is++) {
          SuperkSample_t s = conv::kSKSamples[is];
          if (s == conv::kSKsUndefined) continue;
#ifndef __USE_2D_TEMPLATES__
          exp_data[is] = skexp->CurrPredictedErecoSpectrum(s);
#else
          exp_data[is] = skexp->CurrPredictedErecoVsQ2Spectrum(s);
#endif
          exp_data[is]->SetDirectory(0);
        }//is
        
        //Calculate the chisq for 4 samples
        double chisq = 0;
        
        for (int is = 0 ; is < conv::kns ; is++) {
          SuperkSample_t s = conv::kSKSamples[is];
          if (s == conv::kSKsUndefined || s == conv::kSKsRHCFCFV1RingELike) continue;
          chisq += valor::util::Chisq_Penalty_TH2D(real_data[is],exp_data[is]);
        }
        
        double exponent = -(chisq-scale_toys_based_on_samples_scale)/2;
        double prob = exp(exponent);
        chisq_distribution->Fill(chisq);
        if(runopt->UseMCSampling()){
          double random_number = Rnd::Instance()->RndOther().Uniform(0,1);
          if(random_number <= prob){
            scale_throws_accepted++;
          }
          else{
            scale_throws_rejected++;
            use_it = 0;
          }
        }
        else {
          // Cut out the most unlikely toys, as they will barely contribute to the p-value calculation
          // The p-theta group performed a study (in 2018) and found that cutting out toys with weight < e^{-10} is optimal
          if(exponent < -10){
            scale_throws_rejected++;
            use_it = 0;
          }
          else{
            scale_throws_accepted++;
          }
        }
        if(!use_it) {
          //When the CountFailedToys option is enabled, don't re-throw the toy (even if it was not accepted) because
          //with this option, we're interested in a fixed total number of throws, not a fixed number of accepted toys
          if(!runopt->CountFailedToys()){
            //Redo throw if toy throw fails test
            unphysical_tweak = true;
          }
        }
        else {
          chisq_distribution_accepted->Fill(chisq);
          if(runopt->CountFailedToys()){
            // If CountFailedToys is enabled, we store osc par throws for all toys, but here I need to look
            // at the throws of accepted toys too
            FillOscToyHisto(ph);
          }
        }
        std::cout << "chisq = " << chisq << ", min chisq = " << scale_toys_based_on_samples_scale << ", prob = " << prob << std::endl;
        
        for (int is = 0 ; is < conv::kns ; is++) {
          SuperkSample_t s = conv::kSKSamples[is];
          if (s == conv::kSKsUndefined) continue;
          delete exp_data[is];
        }//is 
      }//Scale_toys_based_on_samples      
    }//applied_negevents_condition
    
    // Fill toy histograms of osc params. Needs to be done after toy is chosen
    util::OscParams::Instance()->FillOscToyHisto(ph);
    util::OscParams::Instance()->Print(cout);
    
    if(!runopt->DoToys()){
      for (int is = 0 ; is < conv::kns ; is++) {
        SuperkSample_t s = conv::kSKSamples[is];
        if (s == conv::kSKsUndefined || s == conv::kSKsRHCFCFV1RingELike) continue;
        
        if(runopt->Using_xml_file()){
          cout << "Use xml file" << endl; // Only read 1st toy
          skobs->ReadSample(xmlfile2.Data(),0,runopt->Using_xml_file(),s);
        }
        else{
          cout << "Use root file" << endl; // Only read 1st toy
          skobs->ReadSample(rootfile2.Data(),0,runopt->Using_xml_file(),s);
        }
      }
    }
    
    //Get the spectra
#ifndef __USE_2D_TEMPLATES__ //1D templates
    TH1D * mock_data[kns];
#else
    TH2D * mock_data[kns];
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
#else
      mock_data[is] = skobs->TwoDSpectrum(s);
#endif
#ifndef __USE_2D_TEMPLATES__
      exp_data[is] = skexp->CurrPredictedErecoSpectrum(s);
#else
      exp_data[is] = skexp->CurrPredictedErecoVsQ2Spectrum(s);
#endif
      exp_data[is]->SetDirectory(0);
      mock_data[is]->SetDirectory(0);
    }//is
    
    // If an XML file was specified, add an entry
    if(runopt->SaveAsXML()) {
      
      // Dataset info
      fxml << "      " << Form("<data_set_%d>",iexperiment) << endl;
      if(runopt->SaveFluctInXML()) {
        fxml << "<!-- ";
        ph.Print(fxml);
        if(runopt->SystTweakMode())
          skexp->PrintSystematicVariations(fxml);
        fxml << " -->" << endl;
      }
      
      // The 'observed' data
      for (int is = 0 ; is < conv::kns ; is++) {
        SuperkSample_t s = conv::kSKSamples[is];
        if (s == conv::kSKsUndefined) continue;
        valor::util::AddEventsToXML(fxml, ftxt, mock_data[is], s, full_points_to_add[is], rand->RndToyMCsyst());
      }//is
      
      fxml << "      " << Form("</data_set_%d>\n",iexperiment) << endl;
      
      fxml << "      " << Form("<exp_data_%d>",iexperiment) << endl;
      
      // The 'observed' data
      for (int is = 0 ; is < conv::kns ; is++) {
        SuperkSample_t s = conv::kSKSamples[is];
        if (s == conv::kSKsUndefined) continue;
        valor::util::AddEventsToXML(fxml, ftxt, exp_data[is], s, full_points_to_add[is], rand->RndToyMCsyst());
      }//is
      
      fxml << "      " << Form("</exp_data_%d>\n",iexperiment) << endl;
      
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
    delete skobs;
  }//iexperiment
  
  //
  // Close files
  //
  if(runopt->SaveAsXML()) {
    fxml << "    </observation>       " << endl;
    fxml << "  </superk>              " << endl;
    fxml << "</analysis_inputs>       ";// << endl;
    fxml.close();
  }
  if(f) {
    util::OscParams::Instance()->Write(f);
    if(runopt->CountFailedToys()){
      WriteAllThrows(f);
    }
    chisq_distribution->Write();
    chisq_distribution_accepted->Write();
    f->Close();
    delete f;
  }//f
  
  cout << endl
  << " - Summary - " << endl
  << "    * Created " << runopt->NExperiments() << " toys" << endl
  << "    * " << weird_unphysical_tweak << " unphysical_tweaks rethrown" << endl
  << "    * " << weird_applied_neg << " applied_negevents_condition rethrown" << endl;
  if(runopt->Scale_toys_based_on_samples())
  {
   cout << "    * " << scale_throws_accepted << " throws accepted" << endl
        << "    * " << scale_throws_rejected << " throws rejected" << endl
    << endl;
  }

  sw.Stop();
  sw.Print();
  delete skexp;
  delete skobs_data;
}


