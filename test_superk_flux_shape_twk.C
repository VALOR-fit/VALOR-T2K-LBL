//
// ------------------------------------------------------------------------------------------------
// Test Cholesky decomposition of SuperK flux covariance matrix
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  const int nexp = 40;

  TH1D * h[nexp];
  TH1D * weight = new TH1D("weights","",100,-2,2);

  for(int iexp=0; iexp<nexp; iexp++) {
    h[iexp] = new TH1D("","",200,0,10);
    skexp.GetBeamInputs()->TweakShape();
    for(int ibin=1; ibin <= h[iexp]->GetNbinsX(); ibin++) {
	double E = h[iexp]->GetBinCenter(ibin);
	double flux = skexp.GetBeamInputs()->SuperkNuMuFlux(E);
        double wght = skexp.GetBeamInputs()->SuperkNuMuFluxShapeWeight(E);
      //cout << "Weight = " << wght << endl;
        weight->Fill(wght);
        h[iexp]->SetBinContent(ibin, flux*wght);
    }
  }

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  for(int iexp=0; iexp < nexp; iexp++) {
    const char * opt = ((iexp==0) ? "" : "same");
    h[iexp]->Draw(opt);
  }
  c->Update();

}

