//
// ------------------------------------------------------------------------------------------------
// Test flux shape tweaking using the Cholesky method.
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

  gSystem->Load("libT2KDisapp2010a");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::numu10a;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  const int nexp = 400;

  TH1D * h[nexp];
  TH1D * weights   = new TH1D("weights",   "flux tweaking weights", 100,0.5,1.5);
  TH2D * weights_e = new TH2D("weights_e", "flux tweaking weights / energy bin",200,0,10, 100, 0.5, 1.5);

  for(int iexp=0; iexp<nexp; iexp++) {
    h[iexp] = new TH1D("","",200,0,10);
    skexp.GetBeamInputs()->TweakShape();
    for(int ibin=1; ibin <= h[iexp]->GetNbinsX(); ibin++) {
	double E = h[iexp]->GetBinCenter(ibin);
	double flux = skexp.GetBeamInputs()->SuperkNuMuFlux(E);
        double wght = skexp.GetBeamInputs()->SuperkNuMuFluxShapeWeight(E);
      //cout << "Weight = " << wght << endl;
        weights->Fill(wght);
        weights_e->Fill(E,wght);
        h[iexp]->SetBinContent(ibin, flux*wght);
    }
  }

  TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  for(int iexp=0; iexp < nexp; iexp++) {
    const char * opt = ((iexp==0) ? "" : "same");
    h[iexp]->Draw(opt);
  }
  c1->Update();

  TCanvas * c2 = new TCanvas("c2","",20,20,500,500);
  c2->Divide(1,2);
  c2->cd(1);
  weights->Draw();
  c2->cd(2);
  weights_e->Draw("colz");
  c2->Update();

}

