
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

void MakeBinning()
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

  // Fill p-theta array
  const double WidthMom = 0.1; // GeV/c
  const double WidthTheta = 10; // degrees
  const int Nmom = 16;
  const int Ntheta = 16;
  double array_mom[Nmom];
  double array_theta[Ntheta];
  
  for(int imom=0;imom<Nmom;imom++){
    array_mom[imom] = (double)imom*(double)WidthMom;
  }
  array_theta[Nmom-1] = 1.5;

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

  //
  
  SuperkMC *skmc = new SuperkMC();

  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined){
      continue;
    }

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
    for(int imom=0;imom<Nmom;imom++){
      for(int itheta=0;itheta<Ntheta;itheta++){
	
	double thlep = array_theta[itheta];
	double thlep_rad = thlep / 180. * TMath::Pi();
	double costhlep = TMath::Cos(thlep_rad);
	
	double plep = array_mom[imom]; // GeV
	double elep = sqrt(plep*plep + mass_lep*mass_lep);
	
	// Get neutrino ereco
	double Ereco = 0;
	if( s==kSKsFCFV1RingELike || s==kSKsFCFV1RingMuLike ){
	  Ereco = skmc->ErecoFunc_Nu(plep*1000.,costhlep,mass_lep*1000.); // Inputs to Ereco are MeV 
	}
	else if( s==kSKsRHCFCFV1RingELike || s==kSKsRHCFCFV1RingMuLike ){
	  Ereco = skmc->ErecoFunc_NuBar(plep*1000.,costhlep,mass_lep*1000.); // Inputs to Ereco are MeV 
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
	double xsecterm = skmc->CalcXSecTerm(Ereco,elep,plep,costhlep,mass_lep,mass_nucl);
	
	// Don't consider {p-theta} combination that are not physical
	if(Ereco < 0.) continue;

	// Fill xsecterm vector
	vec_xsecterm.push_back(xsecterm);

	cout << "p = " << plep 
	     << " / costhlep = " << costhlep 
	     << " --> Ereco = " << Ereco 
	     << " --> XSecTerm = " << xsecterm
	     << endl;
      }
    }
    cout << endl;
    

    // Order the xsecterm vector and make the 1D histogram
    
    sort(vec_xsecterm.begin(), vec_xsecterm.end());

    // Remove identical bin edges (e.g. p=0)

    cout << endl;
    cout << "Vector before erasing identical bin edges:" << endl;
    for(int i=0;i<vec_xsecterm.size();i++){
      cout << vec_xsecterm[i] << " ";
    }
    cout << "Size is " << vec_xsecterm.size() << endl;
    cout << endl;

    // Conditions needed to remove bin edges
    for(int i=0;i<vec_xsecterm.size()-1;i++){   // don't remove the last bin edge 
      if( 
	 ( fabs(vec_xsecterm[i] - vec_xsecterm[i-1]) < 1.E-3 ) ||  // remove if very similar bin edges
	 (vec_xsecterm[i] > 2.)                                    // remove if XSecTerm > 2 (a few events there)
	  ){
	cout << vec_xsecterm[i-1] << " " << vec_xsecterm[i];
	cout << " --> remove element " << i << " --> " << vec_xsecterm[i] << endl;
	vec_xsecterm.erase(vec_xsecterm.begin()+i);
	i=0; // reset and start from the beginning
      }
    }

    cout << endl;
    cout << "Vector after erasing identical bin edges:" << endl;
    for(int i=0;i<vec_xsecterm.size();i++){
      cout << vec_xsecterm[i] << " ";
    }
    cout << "Size is " << vec_xsecterm.size() << endl;
    cout << endl;

    // Get the array
    
    int Nbinedges_xsecterm = vec_xsecterm.size();
    double *array_xsecterm = new double[Nbinedges_xsecterm];
    
    for(int i=0;i<Nbinedges_xsecterm;i++){
      array_xsecterm[i] = vec_xsecterm[i];
      // cout << vec_xsecterm[i] << " " << array_xsecterm[i] << endl;
    }

    
    // Build the histogram
    
    if(s==kSKsFCFV1RingMuLike){
      hXSecTerm_mulike_fhc = new TH1D("hXSecTerm_mulike_fhc","hXSecTerm_mulike_fhc",Nbinedges_xsecterm-1,array_xsecterm);    
      for(int i=0;i<hXSecTerm_mulike_fhc->GetNbinsX();i++){
	hXSecTerm_mulike_fhc->SetBinContent(i+1,hXSecTerm_mulike_fhc->GetNbinsX()-i);
      }
    }
    else if(s==kSKsRHCFCFV1RingMuLike){
      hXSecTerm_mulike_rhc = new TH1D("hXSecTerm_mulike_rhc","hXSecTerm_mulike_rhc",Nbinedges_xsecterm-1,array_xsecterm);    
      for(int i=0;i<hXSecTerm_mulike_rhc->GetNbinsX();i++){
	hXSecTerm_mulike_rhc->SetBinContent(i+1,hXSecTerm_mulike_rhc->GetNbinsX()-i);
      }
    }
    else if(s==kSKsFCFV1RingELike){
      hXSecTerm_elike_fhc = new TH1D("hXSecTerm_elike_fhc","hXSecTerm_elike_fhc",Nbinedges_xsecterm-1,array_xsecterm);    
      for(int i=0;i<hXSecTerm_elike_fhc->GetNbinsX();i++){
	hXSecTerm_elike_fhc->SetBinContent(i+1,hXSecTerm_elike_fhc->GetNbinsX()-i);
      }
    }
    else if(s==kSKsRHCFCFV1RingELike){
      hXSecTerm_elike_rhc = new TH1D("hXSecTerm_elike_rhc","hXSecTerm_elike_rhc",Nbinedges_xsecterm-1,array_xsecterm);    
      for(int i=0;i<hXSecTerm_elike_rhc->GetNbinsX();i++){
	hXSecTerm_elike_rhc->SetBinContent(i+1,hXSecTerm_elike_rhc->GetNbinsX()-i);
      }
    }

  } // Loop over samples
  
  TCanvas *cXSecTerm_mulike_fhc = new TCanvas("cXSecTerm_mulike_fhc","cXSecTerm_mulike_fhc");
  hXSecTerm_mulike_fhc->Draw();
  hXSecTerm_mulike_fhc->GetXaxis()->SetTitle("XSecTerm");
  TCanvas *cXSecTerm_mulike_rhc = new TCanvas("cXSecTerm_mulike_rhc","cXSecTerm_mulike_rhc");
  hXSecTerm_mulike_rhc->Draw();
  hXSecTerm_mulike_rhc->GetXaxis()->SetTitle("XSecTerm");
  TCanvas *cXSecTerm_elike_fhc = new TCanvas("cXSecTerm_elike_fhc","cXSecTerm_elike_fhc");
  hXSecTerm_elike_fhc->Draw();
  hXSecTerm_elike_fhc->GetXaxis()->SetTitle("XSecTerm");
  TCanvas *cXSecTerm_elike_rhc = new TCanvas("cXSecTerm_elike_rhc","cXSecTerm_elike_rhc");
  hXSecTerm_elike_rhc->Draw();
  hXSecTerm_elike_rhc->GetXaxis()->SetTitle("XSecTerm");

  // Print canvases
  cXSecTerm_mulike_fhc->Print("cXSecTerm_mulike_fhc.pdf");
  cXSecTerm_mulike_rhc->Print("cXSecTerm_mulike_rhc.pdf");
  cXSecTerm_elike_fhc->Print("cXSecTerm_elike_fhc.pdf");
  cXSecTerm_elike_rhc->Print("cXSecTerm_elike_rhc.pdf");

  // Print out the binning
  
  cout << endl;
  cout << "Print out the arrays of the bin edges" << endl;
  cout << endl;
  
  int NXSecTerm_mulike_fhc = hXSecTerm_mulike_fhc->GetNbinsX();
  TArrayD *axis_bin_mulike_fhc = (TArrayD*) hXSecTerm_mulike_fhc->GetXaxis()->GetXbins();
  double *array_mulike_fhc = axis_bin_mulike_fhc->GetArray();
  cout << endl;
  cout << "Binning of " << util::AsString(kSKsFCFV1RingMuLike) << endl;
  cout << "NXSecTerm_mulike_fhc = " << NXSecTerm_mulike_fhc << endl;
  cout << "{ ";
  for(int ibin=0;ibin<NXSecTerm_mulike_fhc+1;ibin++){
    cout << array_mulike_fhc[ibin];
    if( ibin<(NXSecTerm_mulike_fhc) ) cout << ", ";
  }
  cout << " }";
  cout << endl;
  
  int NXSecTerm_mulike_rhc = hXSecTerm_mulike_rhc->GetNbinsX();
  TArrayD *axis_bin_mulike_rhc = (TArrayD*) hXSecTerm_mulike_rhc->GetXaxis()->GetXbins();
  double *array_mulike_rhc = axis_bin_mulike_rhc->GetArray();
  cout << endl;
  cout << "Binning of " << util::AsString(kSKsRHCFCFV1RingMuLike) << endl;
  cout << "NXSecTerm_mulike_rhc = " << NXSecTerm_mulike_rhc << endl;
  cout << "{ ";
  for(int ibin=0;ibin<NXSecTerm_mulike_rhc+1;ibin++){
    cout << array_mulike_rhc[ibin];
    if( ibin<(NXSecTerm_mulike_rhc) ) cout << ", ";
  }
  cout << " }";
  cout << endl;

  int NXSecTerm_elike_fhc = hXSecTerm_elike_fhc->GetNbinsX();
  TArrayD *axis_bin_elike_fhc = (TArrayD*) hXSecTerm_elike_fhc->GetXaxis()->GetXbins();
  double *array_elike_fhc = axis_bin_elike_fhc->GetArray();
  cout << endl;
  cout << "Binning of " << util::AsString(kSKsFCFV1RingELike) << endl;
  cout << "NXSecTerm_elike_fhc = " << NXSecTerm_elike_fhc << endl;
  cout << "{ ";
  for(int ibin=0;ibin<NXSecTerm_elike_fhc+1;ibin++){
    cout << array_elike_fhc[ibin];
    if( ibin<(NXSecTerm_elike_fhc) ) cout << ", ";
  }
  cout << " }";
  cout << endl;
  
  int NXSecTerm_elike_rhc = hXSecTerm_elike_rhc->GetNbinsX();
  TArrayD *axis_bin_elike_rhc = (TArrayD*) hXSecTerm_elike_rhc->GetXaxis()->GetXbins();
  double *array_elike_rhc = axis_bin_elike_rhc->GetArray();
  cout << endl;
  cout << "Binning of " << util::AsString(kSKsRHCFCFV1RingELike) << endl;
  cout << "NXSecTerm_elike_rhc = " << NXSecTerm_elike_rhc << endl;
  cout << "{ ";
  for(int ibin=0;ibin<NXSecTerm_elike_rhc+1;ibin++){
    cout << array_elike_rhc[ibin];
    if( ibin<(NXSecTerm_elike_rhc) ) cout << ", ";
  }
  cout << " }";
  cout << endl;

  cout << endl;

  return;
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

