
//
// It takes expected histogram from test_sk_ereco_spectrum_nominal_mainmode_nuebarapp.C output
//

using namespace std;
using namespace valor;
using namespace valor::conv;

void plot_ereco_on_top_of_prediction
(
 const char * sk_sample            = "RHC1Re",                         /* options: "1Rmu", "1Re", "MultiRmu" */
 const char * tag    = "fitmarg_betabar_DEV",    /* tag for output file */
 
 const char * expfile = "sk_ereco_spectrum_RHC1Re_prova-sinsqtheta[12_0.306,13_0.025,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007500,23_0.002509]eV2,betabar_1.000-RHC1Re.root",
 const char * xmlfile = "../../inputs/data_run1to6/skdata_run1to6.xml",
 )
{
  // Get Expectation
  
  TFile *fexp = new TFile(expfile,"READ");
  TH1D *hErecoExp = (TH1D*) fexp->Get("sk_RHC1Retotal");
  
  // Get Data

  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile));
  
  FileStat_t dummy;
  if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) 
    {
      xmlfile2 = xmlfile;
      if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
	cerr<<"Data input file not found: "<<xmlfile<<endl;
	exit(-1);
      }
    }
 
  SuperkObservation * skobs_fakedata = new SuperkObservation;
  bool using_xml_file = true;
  if(using_xml_file){
    skobs_fakedata->Read(xmlfile2.Data(),0,using_xml_file);
  }
  else{
    skobs_fakedata->Read(rootfile2.Data(),0,using_xml_file);
  }
  
  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample


  // Draw 

  TH1D *hErecoData;   
  hErecoData = (TH1D*)skobs_fakedata->ErecoSpectrum(s)->Clone("hErecoData");
  
  hErecoData->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoData->GetYaxis()->SetTitle("Number of events");
  
  hErecoExp->SetLineColor(kBlue);
  
  hErecoData->Draw("pe");
  hErecoExp->Draw("histo same");  

  TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.7);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hErecoExp,"Expected","L");
  leg->AddEntry(hErecoData,"Data","P");
  leg->Draw();
  
}
