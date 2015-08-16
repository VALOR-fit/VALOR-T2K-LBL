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
#include "Conventions/PhysHypothesisId.h"

#include "Build.h"
#include "UtlForMacros.h"
#include "OscParams.h"
#include "RunOpt.h"
#include "Utilities/Env.h"

#include "AnalysisOptions.h"

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

#ifndef __USE_2D_TEMPLATES__ //1D templates
void AddExtraEventsToHist(TH1D * h, string * full_points_to_add);
void AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH1D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand);
#else //2D templates
void AddExtraEventsToHist(TH2D * h, string * full_points_to_add);
void AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH2D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand);
#endif
string * ConstructPointsToAddFromXML(TXMLEngine & xml, XMLNodePointer_t & root_node, valor::conv::SuperkSample_t s, string sextra, int datasetid);

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

int main(int argc, char ** argv)
{
  TStopwatch sw;
  
  Env::Instance()->AutoSet("t2xk");
  cout << *Env::Instance() << endl;
  gROOT->SetBatch(true);
  MyRunOpt runopt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParse3FOscParams);

  TRandom3 rand(runopt.ToySeed()); //random number to see where inside a bin a toy is generated

  TH1D * h_0 = new TH1D ("h_0", "h_0", 100, 0, 10);
  TH1D * h_1 = new TH1D ("h_1", "h_1", 100, 0, 10);
  TH1D * h_2 = new TH1D ("h_2","h_2", 100, 0, 10);
  TH1D * h_3 = new TH1D ("h_3","h_3", 100, 0, 10);
  TH1D * h_4 = new TH1D ("h_4","h_4", 100, 0, 10);
  TH2D * h_4_Escale = new TH2D ("h_4_Escale","h_4_Escale", 100, -5, 5, 100, 0, 10);

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
      full_points_to_add[is] = ConstructPointsToAddFromXML(xml, root_node, s, runopt.ExtraPoints(is), runopt.XMLExtraId());
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

  int weird_applied_neg      = 0;
  int weird_unphysical_tweak = 0;
  for(int iexperiment = runopt.FirstExperiment(), last_experiment = runopt.LastExperiment(); iexperiment <= last_experiment; iexperiment++) {
     cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;

     // Generate mock data
     cout << "- Generating toy-MC..." << endl;
     SuperkObservation * skobs = new SuperkObservation;
     PhysHypothesis ph(runopt.Hypothesis());

     // Set the mass hierarchy
     util::SetMassHierarchy(skexp);

     // Randomize specified osc params for current toy MC expt
     util::OscParams::Instance()->Randomize();

     // Update physics hypothesis with chosen value of oscillation parameters
     util::OscParams::Instance()->UpdatePhysHypothesis(ph);

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
    cerr << endl << "Warning h_0, h_1, ... have not been filled (what are they for?!)" << endl;
    h_0->Write();
    h_1->Write();
    h_2->Write();
    h_3->Write();
    h_4->Write();
    h_4_Escale->Write();
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

string * ConstructPointsToAddFromXML(TXMLEngine & xml, XMLNodePointer_t & root_node, valor::conv::SuperkSample_t s, string sextra, int datasetid)
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
