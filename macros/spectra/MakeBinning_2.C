
//
// ------------------------------------------------------------------------------------------------           
//
// Derive binning of new variables from existing Ereco binning (theta is fixed)
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

#include <vector>       // std::vector

using namespace std;
using namespace valor;
using namespace valor::conv;

double GetXSecTermVal(double SuperkMC *skmc,SuperkSample_t s,double costheta,double momlep /*GeV*/);
double GetCosTheta(double theta_degree);

void MakeBinning_2
(
 const char * sk_sample = "RHC1Re"                         /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
 )
{
  // Define constants
  const double mass_e       = 0.511   / 1000.; // electron mass (GeV)
  const double mass_mu      = 105.658 / 1000.; // muon mass (GeV)
  const double neutronMass  = 939.565 / 1000.; // GeV
  const double protonMass   = 938.272 / 1000.; // GeV
  const double dir_nu[3]    = {0.669764, -0.742179, 0.024223}; // T2K beam direction

  // Define histograms
  TH1D *hXSecTerm_mulike_fhc;
  TH1D *hXSecTerm_mulike_rhc;
  TH1D *hXSecTerm_elike_fhc;
  TH1D *hXSecTerm_elike_rhc;

  TGraph *gBinEdges_Up = new TGraph();
  TGraph *gBinEdges_Low = new TGraph();
  TGraph *gBinEdges_All = new TGraph();

  // Fill p-theta array
  const double WidthMom = 0.05; // 0.1 // GeV/c
  const double WidthTheta = 3; //10; // degrees
  const int Nmom = 31; //16;
  const int Ntheta = 60; //16
  double array_mom[Nmom];
  double array_theta[Ntheta];
  
  for(int imom=0;imom<(Nmom-1);imom++){
    array_mom[imom] = (double)imom*(double)WidthMom;
  }
  array_mom[Nmom-1] = 1.5;

  for(int itheta=0;itheta<(Ntheta-1);itheta++){
    array_theta[itheta] = (double)itheta*(double)WidthTheta;
  }
  array_theta[Ntheta-1] = 180; 

  cout << endl;
  cout << "Momentum bin edges:" << endl;
  for(int imom=0;imom<Nmom;imom++){
    cout << array_mom[imom] << " ";
  }
  cout << endl;
  cout << "Theta bin edges:" << endl;
  for(int itheta=0;itheta<Ntheta;itheta++){
    cout << array_theta[itheta] << " ";
  }
  cout << endl;

  

  double array[22];
  for(int i=0;i<21;i++){
    array[i] = (double)i*(double)0.1;
  }
  array[21] = 5.;
  cout << "XSecTerm bin edges:" << endl;
  cout << "{" << endl;
  for(int i=0;i<22;i++){
    cout << array[i];
    if(i<21) cout << ", ";
  }
  cout << endl;
  cout << "}" << endl;
  cout << endl;

  return;




  
  // Build p-theta histogram
  
  TH2D *hptheta = new TH2D("hptheta","hptheta",Ntheta-1,array_theta,Nmom-1,array_mom);
 
  //
  
  SuperkMC *skmc = new SuperkMC();

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  
  cout << endl;
  cout << util::AsString(s) << endl;
  cout << endl;
  
  vector<double> vec_xsecterm;
  
  // Get lepton mass
  double mass_lep;
  if( s==kSKsRHCFCFV1RingMuLike || s==kSKsFCFV1RingMuLike ){
    mass_lep = mass_mu; // GeV	
  }
  else if( s==kSKsRHCFCFV1RingELike || s==kSKsFCFV1RingELike ){
    mass_lep = mass_e; // GeV	
  }
  else{
    cerr << endl;
    cerr << "SuperkMC::FillTemplateHistogram" << endl;
    cerr << "Using sample " << util::AsString(s) << endl;
    cerr << "Samples other than: " << endl;
    cerr << " - " << util::AsString(kSKsFCFV1RingMuLike) << endl;
    cerr << " - " << util::AsString(kSKsFCFV1RingELike) << endl;
    cerr << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
    cerr << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
    cerr << "are not implemented" << endl;
    exit(1);
  }
  
  // Get nucleon mass (nu-->neutron, nubar-->proton)
  double mass_nucl = 0.;
  if( (s == kSKsRHCFCFV1RingELike) || (s == kSKsRHCFCFV1RingMuLike) ){ // RHC
    mass_nucl = protonMass; // nubar interacts with proton (GeV)
  }
  else if( (s == kSKsFCFV1RingELike) || (s == kSKsFCFV1RingMuLike) ){ // FHC
    mass_nucl = neutronMass; // nubar interacts with neutron (GeV)
  }
  else{
    cerr << endl;
    cerr << "SuperkMC::FillTemplateHistogram" << endl;
    cerr << "Using sample " << util::AsString(s) << endl;
    cerr << "Samples other than: " << endl;
    cerr << " - " << util::AsString(kSKsFCFV1RingMuLike) << endl;
    cerr << " - " << util::AsString(kSKsFCFV1RingELike) << endl;
    cerr << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
    cerr << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
    cerr << "are not implemented" << endl;
    exit(1);
  }
  
  // Fill the XSecTerm binning
  
  int NbinsX = hptheta->GetNbinsX();
  int NbinsY = hptheta->GetNbinsY();
  int NbinsTot = NbinsX * NbinsY;
  
  //for(int ibin=0;ibin<NbinsTot;ibin++){
  for(int imom=0;imom<Nmom;imom++){
    //for(int itheta=10;itheta<13;itheta++){
    for(int itheta=0;itheta<Ntheta;itheta++){
      
      
      // Get neutrino ereco at bin center
      //double plep = array_mom[imom]; // GeV
      //double thlep = array_theta[itheta];      
      //double Ereco = GetErecoVal(skmc,s,costhlep,plep*1000.,mass_lep*1000.); // Inputs to Ereco are MeV 
      //Ereco /= 1000.; // MeV --> GeV
      // if(Ereco < 0.){
      // 	cerr << endl;
      // 	cerr << "Ereco = " << Ereco << endl;
      // 	cerr << "plep = " << plep << endl;
      // 	cerr << "thlep = " << thlep << endl;
      // 	cerr << "costhlep = " << costhlep << endl;
      // 	cerr << "Skip this p-theta values, they are unphysical" << endl;
      // 	cerr << endl;
      // 	continue; //Don't consider {p-theta} combination that are not physical
      // }
      // // Calculate the XSecTerm
      //double xsecterm = GetXSecTermVal(skmc,s,costhlep,plep,mass_lep,mass_nucl); // Inputs are GeV      
      
      
      // Get XSecTerm bin edges 

      double plep = hptheta->GetYaxis()->GetBinCenter(imom+1);
      double thlep = hptheta->GetXaxis()->GetBinCenter(itheta+1);
      double costhlep = GetCosTheta(thlep);
      
      double theta_low    = hptheta->GetXaxis()->GetBinLowEdge(itheta+1);
      double costhlep_low = GetCosTheta(theta_low);
      double plep_low     = hptheta->GetYaxis()->GetBinLowEdge(imom+1);
      double theta_up     = hptheta->GetXaxis()->GetBinUpEdge(itheta+1);
      double costhlep_up  = GetCosTheta(theta_up);
      double plep_up      = hptheta->GetYaxis()->GetBinUpEdge(imom+1);
      
      // low Ereco - low Theta
      
      double XSecTerm_low_low = -999;
      double Ereco_low_low = GetErecoVal(skmc,s,costhlep_low,plep_low*1000.,mass_lep*1000.); //Inputs to Ereco are MeV
      Ereco_low_low /= 1000.; // MeV --> GeV
      if(Ereco_low_low >= 0.) 
	XSecTerm_low_low = GetXSecTermVal(skmc,s,costhlep_low,plep_low,mass_lep,mass_nucl);
      
      double XSecTerm_low_up = -999;
      double Ereco_low_up = GetErecoVal(skmc,s,costhlep_low,plep_up*1000.,mass_lep*1000.); //Inputs to Ereco are MeV
      Ereco_low_up /= 1000.; // MeV --> GeV
      if(Ereco_low_up >= 0.) 
	XSecTerm_low_up = GetXSecTermVal(skmc,s,costhlep_low,plep_up,mass_lep,mass_nucl);
 
      double XSecTerm_up_low = -999;
      double Ereco_up_low = GetErecoVal(skmc,s,costhlep_up,plep_low*1000.,mass_lep*1000.); //Inputs to Ereco are MeV
      Ereco_up_low /= 1000.; // MeV --> GeV
      if(Ereco_up_low >= 0.) 
	XSecTerm_up_low = GetXSecTermVal(skmc,s,costhlep_up,plep_low,mass_lep,mass_nucl);

      double XSecTerm_up_up = -999;
      double Ereco_up_up = GetErecoVal(skmc,s,costhlep_up,plep_up*1000.,mass_lep*1000.); //Inputs to Ereco are MeV
      Ereco_up_up /= 1000.; // MeV --> GeV
      if(Ereco_up_up >= 0.) 
	XSecTerm_up_up = GetXSecTermVal(skmc,s,costhlep_up,plep_up,mass_lep,mass_nucl);
      
      // Get largest bin width
      
      // double XSecTerm_Up = -999;
      // if(Ereco_low_low>0.) XSecTerm_Up = TMath::Max(XSecTerm_Up,XSecTerm_low_low);
      // if(Ereco_low_up>0.)  XSecTerm_Up = TMath::Max(XSecTerm_Up,XSecTerm_low_up);
      // if(Ereco_up_low>0.)  XSecTerm_Up = TMath::Max(XSecTerm_Up,XSecTerm_up_low);
      // if(Ereco_up_up>0.)   XSecTerm_Up = TMath::Max(XSecTerm_Up,XSecTerm_up_up);
      // double XSecTerm_Low = 999;
      // if(Ereco_low_low>0.) XSecTerm_Low = TMath::Min(XSecTerm_Low,XSecTerm_low_low);
      // if(Ereco_low_up>0.)  XSecTerm_Low = TMath::Min(XSecTerm_Low,XSecTerm_low_up);
      // if(Ereco_up_low>0.)  XSecTerm_Low = TMath::Min(XSecTerm_Low,XSecTerm_up_low);
      // if(Ereco_up_up>0.)   XSecTerm_Low = TMath::Min(XSecTerm_Low,XSecTerm_up_up);

      double XSecTerm_Up = -999;
      if(Ereco_low_low>0.) XSecTerm_Up = TMath::Max(XSecTerm_Up,Ereco_low_low);
      if(Ereco_low_up>0.)  XSecTerm_Up = TMath::Max(XSecTerm_Up,Ereco_low_up);
      if(Ereco_up_low>0.)  XSecTerm_Up = TMath::Max(XSecTerm_Up,Ereco_up_low);
      if(Ereco_up_up>0.)   XSecTerm_Up = TMath::Max(XSecTerm_Up,Ereco_up_up);
      double XSecTerm_Low = 999;
      if(Ereco_low_low>0.) XSecTerm_Low = TMath::Min(XSecTerm_Low,Ereco_low_low);
      if(Ereco_low_up>0.)  XSecTerm_Low = TMath::Min(XSecTerm_Low,Ereco_low_up);
      if(Ereco_up_low>0.)  XSecTerm_Low = TMath::Min(XSecTerm_Low,Ereco_up_low);
      if(Ereco_up_up>0.)   XSecTerm_Low = TMath::Min(XSecTerm_Low,Ereco_up_up);
      
      
      int CountBinEdgeNeg = 0;     
      if(Ereco_low_low < 0.) CountBinEdgeNeg++; 
      if(Ereco_low_up < 0.)  CountBinEdgeNeg++;
      if(Ereco_up_low < 0.)  CountBinEdgeNeg++;
      if(Ereco_up_up < 0.)   CountBinEdgeNeg++;
      
      // Skip 
      if(CountBinEdgeNeg>0){
     	cerr << endl;
	cout << "Ereco bin edges: " << Ereco_low_low << " / " << Ereco_low_up << " / " << Ereco_up_low << " / " << Ereco_up_up << endl;
      	cerr << "plep = " << plep << endl;
      	cerr << "thlep = " << thlep << endl;
      	cerr << "costhlep = " << costhlep << endl;
      	cerr << "Skip this p-theta values, they are unphysical" << endl;
      	cerr << endl;
      	continue; //Don't consider {p-theta} combination that are not physical
      }

      if( XSecTerm_Low < 0. && XSecTerm_Up > 0.){
	cerr << endl;
	cerr << "XSecTerm_Low < 0 --> " << XSecTerm_Low << endl;
	cerr << "XSecTerm_Up = " << XSecTerm_Up << endl;
	cerr << "Set XSecTerm_Low to 0" << endl;
      	cerr << endl;
	
	XSecTerm_Low = 0.;

	exit(1);
      }

      cout << endl;
      cout << "Ereco: " << Ereco_low_low << " " << Ereco_low_up << " " << Ereco_up_low << " " << Ereco_up_up 
	   << " / XSecTerm: " << XSecTerm_low_low << " " << XSecTerm_low_up << " " << XSecTerm_up_low << " " << XSecTerm_up_up 
	   << " --> " << XSecTerm_Low << " / " << XSecTerm_Up << endl;

      hptheta->SetBinContent(itheta+1,imom+1,1.); // check if bins are skipped

      
      // Fill TGraph of 2D bin ID Vs XSecTerm bin edge
      
      int BinId2D = hptheta->GetBin(itheta+1,imom+1);
      
      gBinEdges_Up->SetPoint (gBinEdges_Up->GetN(), XSecTerm_Up, BinId2D);
      gBinEdges_Low->SetPoint(gBinEdges_Low->GetN(),XSecTerm_Low,BinId2D);
      
      gBinEdges_All->SetPoint(gBinEdges_All->GetN(),XSecTerm_Low,BinId2D);
      gBinEdges_All->SetPoint(gBinEdges_All->GetN(),XSecTerm_Up, BinId2D);

      continue;


      // Fill xsecterm vector
      vec_xsecterm.push_back(xsecterm);
      
      cout << "p = " << plep 
	   << " / costhlep = " << costhlep 
	   << " --> Ereco = " << Ereco
	   << " --> XSecTerm = " << xsecterm
	   << " --> XSecTerm my = " << GetXSecTermVal(skmc,s,costhlep,plep,mass_lep,mass_nucl)
	   << endl;
      
    } // loop over theta bins
  } // loop over momentum bins
  cout << endl;
  
  TCanvas *cptheta = new TCanvas("cptheta","cptheta");
  hptheta->GetXaxis()->SetTitle("#theta (degrees)");
  hptheta->GetYaxis()->SetTitle("p (GeV)");
  hptheta->Draw("colz");
  
  TCanvas *cBinEdges_All = new TCanvas("cBinEdges_All","cBinEdges_All");
  gBinEdges_All->GetXaxis()->SetTitle("E_{Reco}");
  gBinEdges_All->GetYaxis()->SetTitle("{p,#theta} bin ID");
  gBinEdges_Low->SetMarkerColor(kBlue);
  //gBinEdges_Up ->SetMarkerColor(KRed);
  gBinEdges_All->Draw("ap");
  gBinEdges_Low->Draw("p same");
  gBinEdges_Up->Draw("p same");



  
  // // Order the xsecterm vector and make the 1D histogram
  
  // sort(vec_xsecterm.begin(), vec_xsecterm.end());
  
  // // Remove identical bin edges (e.g. p=0)
  
  // cout << endl;
  // cout << "Vector before erasing identical bin edges:" << endl;
  // for(int i=0;i<vec_xsecterm.size();i++){
  //   cout << vec_xsecterm[i] << " ";
  // }
  // cout << "Size is " << vec_xsecterm.size() << endl;
  // cout << endl;
  
  // // Conditions needed to remove bin edges
  // for(int i=0;i<vec_xsecterm.size()-1;i++){   // don't remove the last bin edge 
  //   if( 
  //      ( fabs(vec_xsecterm[i] - vec_xsecterm[i-1]) < 1.E-3 ) ||  // remove if very similar bin edges
  //      (vec_xsecterm[i] > 2.)                                    // remove if XSecTerm > 2 (a few events there)
  // 	){
  //     cout << vec_xsecterm[i-1] << " " << vec_xsecterm[i];
  //     cout << " --> remove element " << i << " --> " << vec_xsecterm[i] << endl;
  //     vec_xsecterm.erase(vec_xsecterm.begin()+i);
  //     i=0; // reset and start from the beginning
  //   }
  // }
  
  // cout << endl;
  // cout << "Vector after erasing identical bin edges:" << endl;
  // for(int i=0;i<vec_xsecterm.size();i++){
  //   cout << vec_xsecterm[i] << " ";
  // }
  // cout << "Size is " << vec_xsecterm.size() << endl;
  // cout << endl;
  
  // // Get the array
  
  // int Nbinedges_xsecterm = vec_xsecterm.size();
  // double *array_xsecterm = new double[Nbinedges_xsecterm];
  
  // for(int i=0;i<Nbinedges_xsecterm;i++){
  //   array_xsecterm[i] = vec_xsecterm[i];
  //   // cout << vec_xsecterm[i] << " " << array_xsecterm[i] << endl;
  // }
  
  
  // // Build the histogram
  
  // if(s==kSKsFCFV1RingMuLike){
  //   hXSecTerm_mulike_fhc = new TH1D("hXSecTerm_mulike_fhc","hXSecTerm_mulike_fhc",Nbinedges_xsecterm-1,array_xsecterm);    
  //   for(int i=0;i<hXSecTerm_mulike_fhc->GetNbinsX();i++){
  //     hXSecTerm_mulike_fhc->SetBinContent(i+1,hXSecTerm_mulike_fhc->GetNbinsX()-i);
  //   }
  // }
  // else if(s==kSKsRHCFCFV1RingMuLike){
  //   hXSecTerm_mulike_rhc = new TH1D("hXSecTerm_mulike_rhc","hXSecTerm_mulike_rhc",Nbinedges_xsecterm-1,array_xsecterm);    
  //   for(int i=0;i<hXSecTerm_mulike_rhc->GetNbinsX();i++){
  //     hXSecTerm_mulike_rhc->SetBinContent(i+1,hXSecTerm_mulike_rhc->GetNbinsX()-i);
  //   }
  // }
  // else if(s==kSKsFCFV1RingELike){
  //   hXSecTerm_elike_fhc = new TH1D("hXSecTerm_elike_fhc","hXSecTerm_elike_fhc",Nbinedges_xsecterm-1,array_xsecterm);    
  //   for(int i=0;i<hXSecTerm_elike_fhc->GetNbinsX();i++){
  //     hXSecTerm_elike_fhc->SetBinContent(i+1,hXSecTerm_elike_fhc->GetNbinsX()-i);
  //   }
  // }
  // else if(s==kSKsRHCFCFV1RingELike){
  //   hXSecTerm_elike_rhc = new TH1D("hXSecTerm_elike_rhc","hXSecTerm_elike_rhc",Nbinedges_xsecterm-1,array_xsecterm);    
  //   for(int i=0;i<hXSecTerm_elike_rhc->GetNbinsX();i++){
  //     hXSecTerm_elike_rhc->SetBinContent(i+1,hXSecTerm_elike_rhc->GetNbinsX()-i);
  //   }
  // }
  
  // //  } // Loop over samples
  



  // TCanvas *cXSecTerm_mulike_fhc = new TCanvas("cXSecTerm_mulike_fhc","cXSecTerm_mulike_fhc");
  // hXSecTerm_mulike_fhc->Draw();
  // hXSecTerm_mulike_fhc->GetXaxis()->SetTitle("XSecTerm");
  // TCanvas *cXSecTerm_mulike_rhc = new TCanvas("cXSecTerm_mulike_rhc","cXSecTerm_mulike_rhc");
  // hXSecTerm_mulike_rhc->Draw();
  // hXSecTerm_mulike_rhc->GetXaxis()->SetTitle("XSecTerm");
  // TCanvas *cXSecTerm_elike_fhc = new TCanvas("cXSecTerm_elike_fhc","cXSecTerm_elike_fhc");
  // hXSecTerm_elike_fhc->Draw();
  // hXSecTerm_elike_fhc->GetXaxis()->SetTitle("XSecTerm");
  // TCanvas *cXSecTerm_elike_rhc = new TCanvas("cXSecTerm_elike_rhc","cXSecTerm_elike_rhc");
  // hXSecTerm_elike_rhc->Draw();
  // hXSecTerm_elike_rhc->GetXaxis()->SetTitle("XSecTerm");

  // // Print canvases
  // cXSecTerm_mulike_fhc->Print("cXSecTerm_mulike_fhc.pdf");
  // cXSecTerm_mulike_rhc->Print("cXSecTerm_mulike_rhc.pdf");
  // cXSecTerm_elike_fhc->Print("cXSecTerm_elike_fhc.pdf");
  // cXSecTerm_elike_rhc->Print("cXSecTerm_elike_rhc.pdf");

  // // Print out the binning
  
  // cout << endl;
  // cout << "Print out the arrays of the bin edges" << endl;
  // cout << endl;
  
  // int NXSecTerm_mulike_fhc = hXSecTerm_mulike_fhc->GetNbinsX();
  // TArrayD *axis_bin_mulike_fhc = (TArrayD*) hXSecTerm_mulike_fhc->GetXaxis()->GetXbins();
  // double *array_mulike_fhc = axis_bin_mulike_fhc->GetArray();
  // cout << endl;
  // cout << "Binning of " << util::AsString(kSKsFCFV1RingMuLike) << endl;
  // cout << "NXSecTerm_mulike_fhc = " << NXSecTerm_mulike_fhc << endl;
  // cout << "{ ";
  // for(int ibin=0;ibin<NXSecTerm_mulike_fhc+1;ibin++){
  //   cout << array_mulike_fhc[ibin];
  //   if( ibin<(NXSecTerm_mulike_fhc) ) cout << ", ";
  // }
  // cout << " }";
  // cout << endl;
  
  // int NXSecTerm_mulike_rhc = hXSecTerm_mulike_rhc->GetNbinsX();
  // TArrayD *axis_bin_mulike_rhc = (TArrayD*) hXSecTerm_mulike_rhc->GetXaxis()->GetXbins();
  // double *array_mulike_rhc = axis_bin_mulike_rhc->GetArray();
  // cout << endl;
  // cout << "Binning of " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
  // cout << "NXSecTerm_mulike_rhc = " << NXSecTerm_mulike_rhc << endl;
  // cout << "{ ";
  // for(int ibin=0;ibin<NXSecTerm_mulike_rhc+1;ibin++){
  //   cout << array_mulike_rhc[ibin];
  //   if( ibin<(NXSecTerm_mulike_rhc) ) cout << ", ";
  // }
  // cout << " }";
  // cout << endl;

  // int NXSecTerm_elike_fhc = hXSecTerm_elike_fhc->GetNbinsX();
  // TArrayD *axis_bin_elike_fhc = (TArrayD*) hXSecTerm_elike_fhc->GetXaxis()->GetXbins();
  // double *array_elike_fhc = axis_bin_elike_fhc->GetArray();
  // cout << endl;
  // cout << "Binning of " << util::AsString(kSKsFCFV1RingELike) << endl;
  // cout << "NXSecTerm_elike_fhc = " << NXSecTerm_elike_fhc << endl;
  // cout << "{ ";
  // for(int ibin=0;ibin<NXSecTerm_elike_fhc+1;ibin++){
  //   cout << array_elike_fhc[ibin];
  //   if( ibin<(NXSecTerm_elike_fhc) ) cout << ", ";
  // }
  // cout << " }";
  // cout << endl;
  
  // int NXSecTerm_elike_rhc = hXSecTerm_elike_rhc->GetNbinsX();
  // TArrayD *axis_bin_elike_rhc = (TArrayD*) hXSecTerm_elike_rhc->GetXaxis()->GetXbins();
  // double *array_elike_rhc = axis_bin_elike_rhc->GetArray();
  // cout << endl;
  // cout << "Binning of " << util::AsString(kSKsRHCFCFV1RingELike) << endl;
  // cout << "NXSecTerm_elike_rhc = " << NXSecTerm_elike_rhc << endl;
  // cout << "{ ";
  // for(int ibin=0;ibin<NXSecTerm_elike_rhc+1;ibin++){
  //   cout << array_elike_rhc[ibin];
  //   if( ibin<(NXSecTerm_elike_rhc) ) cout << ", ";
  // }
  // cout << " }";
  // cout << endl;

  // cout << endl;

  // return;
}






//double Q2 = this->CalcQ2(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep); // Q^2   
//double Q2 = this->CalcQd2(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep); // Qd^2   
//double Q2 = this->CalcEhad(ereco_nu,emax_lep); 
//double Q2 = this->CalcInelasticity(ereco_nu,emax_lep);
//double Q2 = this->CalcBjorkenScale(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep,mass_nucl); 
//double Q2 = this->CalcWhadr(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep,mass_nucl);  
//double Q2 = this->CalcXSecTerm(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep,mass_nucl);  

//const int nbins_ereco = util::NumErecoBins(s);
//const double * ereco_bin_edges = util::ErecoBinEdges(s);
//cout << util::AsString(s) << endl;
//cout << "# of Ereco bins = " << nbins_ereco << endl;
//cout << "Ereco bin edges:" << endl;
//for(int ibin=0;ibin<nbins_ereco;ibin++){
// Get lepton momentum;
//double Mom_lep  = ereco_bin_edges[ibin] * 1000.; // MeV





double GetXSecTermVal( SuperkMC *skmc,SuperkSample_t s,
		       double costheta,double momlep /*GeV*/,
		       double mass_lep /*GeV*/, double mass_nucl /*GeV*/
		       ){
  
  double elep = sqrt(momlep*momlep + mass_lep*mass_lep); // 
  
  // Get neutrino ereco
  double Ereco = 0;
  if( s==kSKsFCFV1RingELike || s==kSKsFCFV1RingMuLike ){
    Ereco = skmc->ErecoFunc_Nu(momlep*1000.,costheta,mass_lep*1000.); // Inputs to Ereco are MeV 
  }
  else if( s==kSKsRHCFCFV1RingELike || s==kSKsRHCFCFV1RingMuLike ){
    Ereco = skmc->ErecoFunc_NuBar(momlep*1000.,costheta,mass_lep*1000.); // Inputs to Ereco are MeV 
  }
  else{
    cerr << endl;
    cerr << "SuperkMC::FillTemplateHistogram" << endl;
    cerr << "Using sample " << util::AsString(s) << endl;
    cerr << "Samples other than: " << endl;
    cerr << " - " << util::AsString(kSKsFCFV1RingMuLike) << endl;
    cerr << " - " << util::AsString(kSKsFCFV1RingELike) << endl;
    cerr << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
    cerr << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
    cerr << "are not implemented" << endl;
    exit(1);
  }
  Ereco /= 1000.; // MeV --> GeV

  // Calculate the XSecTerm
  double xsecterm = skmc->CalcXSecTerm(Ereco,elep,momlep,costheta,mass_lep,mass_nucl);

  // Don't consider {p-theta} combination that are not physical
  // if(Ereco < 0.){
  //   cerr << endl;
  //   cerr << "Ereco is negative: " << Ereco << endl;
  //   cerr << endl;
  //   //return -999;
  // }
  //else           return xsecterm;

  return xsecterm;
}


double GetErecoVal( SuperkMC *skmc, SuperkSample_t s,
		    double costheta, double momlep /*MeV*/,
		    double mass_lep /*MeV*/
		    ){
  double Ereco = 0;
  if( s==kSKsFCFV1RingELike || s==kSKsFCFV1RingMuLike ){
    Ereco = skmc->ErecoFunc_Nu(momlep,costheta,mass_lep); 
  }
  else if( s==kSKsRHCFCFV1RingELike || s==kSKsRHCFCFV1RingMuLike ){
    Ereco = skmc->ErecoFunc_NuBar(momlep,costheta,mass_lep); 
  }
  else{
    cerr << endl;
    cerr << "SuperkMC::FillTemplateHistogram" << endl;
    cerr << "Using sample " << util::AsString(s) << endl;
    cerr << "Samples other than: " << endl;
    cerr << " - " << util::AsString(kSKsFCFV1RingMuLike) << endl;
    cerr << " - " << util::AsString(kSKsFCFV1RingELike) << endl;
    cerr << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
    cerr << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
    cerr << "are not implemented" << endl;
    exit(1);
  }
  return Ereco;
}

double GetCosTheta(double theta_degree){
  double theta_rad = theta_degree / 180. * TMath::Pi();
  double costheta = TMath::Cos(theta_rad);
  return costheta;
}
