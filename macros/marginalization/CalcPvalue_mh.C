
//
// ------------------------------------------------------------------------------------------------           
//
// Read and analyze the output from run_mh_marginalization.cxx
// Compute p-values using deltachisq01_glob
// Store final results in a output root file to be used by other macros
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

#include <vector>

using namespace std;
using namespace valor;
using namespace valor::conv;

void CalcPvalue_mh
(
 // Fit of Run 1-6. Best-fit = Normal Hierarchy
 //TString mytag = "Pvalue_IH_10kToys_Run16",
 TString mytag = "Pvalue_NH_10kToys_Asimov",

 //double deltachisq_nh_ih_data = -2.0816342, // chisq(NH)-chisq(IH) - Run 1-6 data set 
 double deltachisq_nh_ih_data = -1.0468992, // Asimov A (Normal Hierarchy)
 
 char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_nh_normal_nexp_10k_beta1.root" // 10k fake data sets, 100k marg toys, Normal Hierarchy
 //char infilename[300] = "../TECH_NOTE_14_3_16/MASSHIERARCHY/FIT/fit_fitmarg_masshier_10k_ih_inverted_nexp_10k_beta1.root" // 10k fake data sets, 100k marg toys, Inverted Hierarchy
 )
{

  // Set data Feldman-Cousins deltachisq later from deltachisq_nh_ih 
  double deltachisq_fc_data = conv::kBadResult; // chisq(True)-chisq(BF) - Run 1-6 data set 
  
  // const double MaxHisto = 20.;
  // const double MinHisto = -20.;
  // const double NBinHisto = 200.;

  const double MaxHisto = 20.;
  const double MinHisto = -20.;
  const double NBinHisto = 200.;
  
  char name[200];

  TFile *filein = new TFile(infilename,"READ");
  TTree *treein = (TTree*)filein->Get("fitnt");
  
  vector<double> vec_chisq_nh_ih;
  vector<double> vec_chisq_fc;
  
  // Tree variables. Set globals because used by both trees
  bool IsFit2D_nh=false;
  bool IsFit2D_ih=false;
  int run=-999;
  string *tag;
  int nexpts=-999;
  int ToyID=-999;
  double truepar0=-999;
  double truepar1=-999;    
  string *trueMH;
  string *trueparname0;
  string *trueparname1;
  string *MH_glob;
  double deltachisq_fc_glob=conv::kALargeNum;  
  double deltachisq_nh_ih_glob=conv::kALargeNum;  
  double chisqmin_glob=conv::kALargeNum;
  double chisqmin_nh=conv::kALargeNum;
  double chisqmin_ih=conv::kALargeNum;
  TH2D *hErecoToy_1Rmu = new TH2D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_1Re = new TH2D("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Rmu = new TH2D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Re = new TH2D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_sample[conv::kns] = {};
  TH2D *hErecoExp_1Rmu = new TH2D("hErecoExp2D_1Rmu","hErecoExp2D_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_1Re = new TH2D("hErecoExp2D_1Re","hErecoExp2D_1Re",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Rmu = new TH2D("hErecoExp2D_RHC1Rmu","hErecoExp2D_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Re = new TH2D("hErecoExp2D_RHC1Re","hErecoExp2D_RHC1Re",1,0,1,1,0,1); 

  treein->SetBranchAddress("ToyID",&ToyID);
  treein->SetBranchAddress("Job",&run);
  treein->SetBranchAddress("Tag",&tag);
  treein->SetBranchAddress("nexpts",&nexpts);
  treein->SetBranchAddress("truepar0", &truepar0);
  treein->SetBranchAddress("truepar1", &truepar1);
  treein->SetBranchAddress("trueMH",   &trueMH);
  treein->SetBranchAddress("trueparname0",&trueparname0);
  treein->SetBranchAddress("trueparname1",&trueparname1);
  treein->SetBranchAddress("deltachisq0BF_glob",&deltachisq_fc_glob);
  treein->SetBranchAddress("deltachisq01_glob",&deltachisq_nh_ih_glob);
  treein->SetBranchAddress("MH_glob",   &MH_glob);
  treein->SetBranchAddress("chisqmin_glob",&chisqmin_glob);
  treein->SetBranchAddress("chisqmin_nh",&chisqmin_nh);
  treein->SetBranchAddress("chisqmin_ih",&chisqmin_ih);
  treein->SetBranchAddress("hErecoToy_1Rmu",     &hErecoToy_1Rmu     );
  treein->SetBranchAddress("hErecoToy_1Re",      &hErecoToy_1Re      );
  treein->SetBranchAddress("hErecoToy_RHC1Rmu",  &hErecoToy_RHC1Rmu  );
  treein->SetBranchAddress("hErecoToy_RHC1Re",   &hErecoToy_RHC1Re   );
  treein->SetBranchAddress("hErecoExp_1Rmu",     &hErecoExp_1Rmu     );
  treein->SetBranchAddress("hErecoExp_1Re",      &hErecoExp_1Re      );
  treein->SetBranchAddress("hErecoExp_RHC1Rmu",  &hErecoExp_RHC1Rmu  );
  treein->SetBranchAddress("hErecoExp_RHC1Re",   &hErecoExp_RHC1Re   );
  
  // Declare histograms
  TH1D *hnevents_rhc1re = new TH1D("hnevents_rhc1re","hnevents_rhc1re",NBinHisto,0.,MaxHisto);  
  TH1D *hnevents_1re = new TH1D("hnevents_1re","hnevents_1re",NBinHisto,0.,MaxHisto);  
  TH1D *hnevents_rhc1rmu = new TH1D("hnevents_rhc1rmu","hnevents_rhc1rmu",NBinHisto,0.,MaxHisto);  
  TH1D *hnevents_1rmu = new TH1D("hnevents_1rmu","hnevents_1rmu",NBinHisto,0.,MaxHisto);  
  
  TH1D *hdeltachisq_nh_ih = new TH1D("hdeltachisq_nh_ih","hdeltachisq_nh_ih",NBinHisto,MinHisto,MaxHisto);
  TH1D *hdeltachisq_fc = new TH1D("hdeltachisq_fc","hdeltachisq_fc",NBinHisto,MinHisto,MaxHisto);
  TH1D *hchisq_glob = new TH1D("hchisq_glob","hchisq_glob",400,0,400);
  TH1D *hchisq_nh = new TH1D("hchisq_nh","hchisq_nh",400,0,400);
  TH1D *hchisq_ih = new TH1D("hchisq_ih","hchisq_ih",400,0,400);
  TH1D *hMH_glob = new TH1D("hMH_glob","hMH_glob",2,0.,2.);
  
  // loop over the tree  
  
  int Nentries = treein->GetEntries();
  for(int ientry=0;ientry<Nentries;ientry++){
  
    treein->GetEntry(ientry);

    if(!(ientry%500)) cout << "Entry " << ientry << endl;

    double nevents_rhc1re = (double)hErecoToy_RHC1Re->Integral();    
    double nevents_1re = (double)hErecoToy_1Re->Integral();    
    double nevents_rhc1rmu = (double)hErecoToy_RHC1Rmu->Integral();    
    double nevents_1rmu = (double)hErecoToy_1Rmu->Integral();    
    if(
       (nevents_rhc1re<0) ||
       (nevents_1re<0) ||
       (nevents_rhc1rmu<0) ||
       (nevents_1rmu<0) 
       ){
      cerr << endl;
      cerr << "ToyID: " << ToyID << endl;
      cerr << "# of events RHC1Re  = " << nevents_rhc1re  << " < 0" << endl;
      cerr << "# of events 1Re     = " << nevents_1re     << " < 0" << endl;
      cerr << "# of events RHC1Rmu = " << nevents_rhc1rmu << " < 0" << endl;
      cerr << "# of events 1Rmu    = " << nevents_1rmu    << " < 0" << endl;
      cerr << endl;
      exit(1);
    }

    hnevents_rhc1re->Fill(nevents_rhc1re);    
    hnevents_1re->Fill(nevents_1re);    
    hnevents_rhc1rmu->Fill(nevents_rhc1rmu);    
    hnevents_1rmu->Fill(nevents_1rmu);
    
    // cout << "chisqmin_nh = " << chisqmin_nh
    //  	 << " | chisqmin_ih = " << chisqmin_ih
    //  	 << " | chisqmin_glob = " << chisqmin_glob
    // 	 << " | deltachisq_nh_ih_glob = " << deltachisq_nh_ih_glob 
    // 	 << " | deltachisq_fc_glob = " << deltachisq_fc_glob 
    //  	 << " |  bestfit MH = " << *MH_glob 
    // 	 << endl;
    
    double my_deltachisq_nh_ih = chisqmin_nh - chisqmin_ih;
    if(my_deltachisq_nh_ih != deltachisq_nh_ih_glob){
      cout << "my_deltachisq_nh_ih = " << my_deltachisq_nh_ih << endl;
      cout << "deltachisq_nh_ih_glob = " << deltachisq_nh_ih_glob << endl;
      exit(1);
    }
    hdeltachisq_nh_ih->Fill(deltachisq_nh_ih_glob);
    hdeltachisq_fc->Fill(deltachisq_fc_glob);

    // Minimized chisq distribution
    hchisq_glob->Fill(chisqmin_glob);
    hchisq_nh->Fill(chisqmin_nh);
    hchisq_ih->Fill(chisqmin_ih);
          
    if(*MH_glob=="normal"){
      hMH_glob->Fill(0.1);
    }
    else if(*MH_glob=="inverted"){
      hMH_glob->Fill(1.1);
    }
    else{
      cerr << endl;
      cerr << "MH_glob = " << *MH_glob << endl;
      cerr << "Best-fit Mass Hierarchy has an incorrect value!" << endl;
      cerr << endl;
      exit(1);
    }

    // Fill vectors for pvalue calculations
    vec_chisq_nh_ih.push_back(deltachisq_nh_ih_glob);
    vec_chisq_fc.push_back(deltachisq_fc_glob);

  } // end loop over toys


  if(deltachisq_nh_ih_glob<-5 || deltachisq_nh_ih_glob>20){
    cout << "deltachisq_nh_ih_glob = " << deltachisq_nh_ih_glob << endl;
    exit(1);
  }
  if(deltachisq_fc_glob<0.){
    cout << "deltachisq_fc_glob = " << deltachisq_fc_glob << endl;
    exit(1);
  }

  
  //
  // Check the Data deltachisq value, depending on the True Mass Hierarchy in the toys
  // If trueMH = MH_glob --> set deltachisq_fc_data to 0
  // Otherwise --> set deltachisq_fc_data to |deltachisq_nh_ih_data|
  //
  
  string databfMH = ""; // best-fit Mass Hierarchy from fit - Run 1-6 data
  if(deltachisq_nh_ih_data < 0) databfMH = "normal";

  cout << endl;
  cout << "//////////// WARNING ////////////" << endl;
  if(*trueMH == databfMH){
    cout << endl;
    cout << "In Data set:" << endl;
    cout << "True Mass Hierarchy = Best-Fit Data Mass Hierarchy" << endl;
    cout << "Set deltachisq_fc_data (Feldman-Cousins) to 0" << endl;
    cout << endl;    
    deltachisq_fc_data = 0.;
  }
  else{
    cout << endl;
    cout << "In Data set:" << endl;
    cout << "True Mass Hierarchy != Best-Fit Data Mass Hierarchy" << endl;
    cout << "Set deltachisq_fc_data (Feldman-Cousins) to |chisq(NH)-chisq(IH)|" << endl;
    cout << endl;
    deltachisq_fc_data = fabs(deltachisq_nh_ih_data);
  }
  cout << "////////////////////////////////" << endl;
  cout << endl;
  
  

  // Compute p-values (NH-IH)

  double Ntoys = (double)Nentries;
  double pvalue_nh_ih_right = GetPvalue_rightside(vec_chisq_nh_ih,Ntoys,deltachisq_nh_ih_data);
  double pvalue_nh_ih_left = GetPvalue_leftside(vec_chisq_nh_ih,Ntoys,deltachisq_nh_ih_data);
  double pvalue_nh_ih_glob = -999;
  
  if(*trueMH == "normal"){
    pvalue_nh_ih_glob = pvalue_nh_ih_right;
  }
  else if(*trueMH == "inverted"){    
    pvalue_nh_ih_glob = pvalue_nh_ih_left;    
  }
  else{
    cout << endl;
    cout << "The true Mass Hierarchy is not known: " << *trueMH << endl;
    cout << endl;
  }

  // Calculate the quantile (p-value = gauss area of both tails --> divide by 2)
  TF1 f2("f2","ROOT::Math::gaussian_quantile(1.-(x/2.),1.) ",0.,1.);
  double zscore_nh_ih_glob = f2.Eval(pvalue_nh_ih_glob); 

  // Compute p-values (a-la-Feldman-Cousins)
  double pvalue_fc_glob = GetPvalue_rightside(vec_chisq_fc,Ntoys,deltachisq_fc_data);
  double zscore_fc_glob = f2.Eval(pvalue_fc_glob); 
  


  
  // Draw histograms
  
  // treein->GetEntry(493);
  // //treein->GetEntry(4739);
  // TCanvas *cErecoToy_GetEntry = new TCanvas("cErecoToy_GetEntry","cErecoToy_GetEntry");
  // hErecoToy_RHC1Re->DrawCopy();

  TCanvas *cMH = new TCanvas("cMH","cMH");
  hMH_glob->GetYaxis()->SetTitle("Number of entries");
  hMH_glob->GetXaxis()->SetTitle("Best-fit Mass Hierarchy");
  hMH_glob->DrawCopy();
  hMH_glob->Draw("same");

  TCanvas *cchisq_nh = new TCanvas("cchisq_nh","cchisq_nh");
  hchisq_nh->GetYaxis()->SetTitle("Number of entries");
  hchisq_nh->GetXaxis()->SetTitle("#chi^{2} (Normal Hierarchy)");
  hchisq_nh->DrawCopy();
  hchisq_nh->Draw("same");
  //cchisq_nh->SetLogy();

  TCanvas *cchisq_ih = new TCanvas("cchisq_ih","cchisq_ih");
  hchisq_ih->GetYaxis()->SetTitle("Number of entries");
  hchisq_ih->GetXaxis()->SetTitle("#chi^{2} (Inverted Hierarchy)");
  hchisq_ih->DrawCopy();
  hchisq_ih->Draw("same");
  //cchisq_ih->SetLogy();

  TCanvas *cchisq_glob = new TCanvas("cchisq_glob","cchisq_glob");
  hchisq_glob->GetYaxis()->SetTitle("Number of entries");
  hchisq_glob->GetXaxis()->SetTitle("#chi^{2}_{glob}");
  hchisq_glob->DrawCopy();
  hchisq_glob->Draw("same");
  //cchisq_glob->SetLogy();

  TCanvas *cPvalue_NH_IH = new TCanvas("cPvalue_NH_IH","cPvalue_NH_IH");
  //cPvalue->SetLogy();
  hdeltachisq_nh_ih->GetYaxis()->SetTitle("Number of entries");
  hdeltachisq_nh_ih->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}(NH) - #chi^{2}(IH)");
  hdeltachisq_nh_ih->DrawCopy();  
  TLine *ldeltachisq_nh_ih = new TLine(deltachisq_nh_ih_data,0.,deltachisq_nh_ih_data,hdeltachisq_nh_ih->GetMaximum());
  ldeltachisq_nh_ih->SetLineColor(kRed);
  ldeltachisq_nh_ih->SetLineWidth(2);
  ldeltachisq_nh_ih->Draw();
  TLatex *latPvalue_nh_ih = new TLatex(hdeltachisq_nh_ih->GetXaxis()->GetBinLowEdge(5),hdeltachisq_nh_ih->GetMaximum()*0.7, TString::Format("P-value = %1.3f",pvalue_nh_ih_glob).Data() );
  latPvalue_nh_ih->Draw("same");
  latPvalue_nh_ih->SetTextSize(0.06);
  TLine *lPvalue_nh_ih = new TLine(pvalue_nh_ih_glob,pvalue_nh_ih_glob,pvalue_nh_ih_glob,pvalue_nh_ih_glob);  
  TLatex *latZscore_nh_ih = new TLatex(hdeltachisq_nh_ih->GetXaxis()->GetBinLowEdge(5),hdeltachisq_nh_ih->GetMaximum()*0.5, TString::Format("Z = %1.2f#sigma",zscore_nh_ih_glob).Data() );
  latZscore_nh_ih->Draw("same");
  latZscore_nh_ih->SetTextSize(0.06);
  TLine *lZscore_nh_ih = new TLine(zscore_nh_ih_glob,zscore_nh_ih_glob,zscore_nh_ih_glob,zscore_nh_ih_glob);
  TLegend *legPvalue_nh_ih = new TLegend(0.5,0.6,0.8,0.8);
  legPvalue_nh_ih->SetFillColor(0);
  legPvalue_nh_ih->SetLineColor(0);
  legPvalue_nh_ih->SetTextFont(132);
  legPvalue_nh_ih->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_nh_ih_data); 
  legPvalue_nh_ih->AddEntry(ldeltachisq_nh_ih,name,"l");
  legPvalue_nh_ih->Draw();  

  TCanvas *cPvalue_FC = new TCanvas("cPvalue_FC","cPvalue_FC");
  //cPvalue->SetLogy();
  hdeltachisq_fc->GetYaxis()->SetTitle("Number of entries");
  hdeltachisq_fc->GetXaxis()->SetTitle("#Delta #chi^{2} = #chi^{2}_{True} - #chi^{2}_{Best-Fit}");
  hdeltachisq_fc->DrawCopy();  
  TLine *ldeltachisq_fc = new TLine(deltachisq_fc_data,0.,deltachisq_fc_data,hdeltachisq_fc->GetMaximum());
  ldeltachisq_fc->SetLineColor(kRed);
  ldeltachisq_fc->SetLineWidth(2);
  ldeltachisq_fc->Draw();
  TLatex *latPvalue_fc = new TLatex(hdeltachisq_fc->GetXaxis()->GetBinLowEdge(5),hdeltachisq_fc->GetMaximum()*0.7, TString::Format("P-value = %1.3f",pvalue_fc_glob).Data() );
  latPvalue_fc->Draw("same");
  latPvalue_fc->SetTextSize(0.06);
  TLine *lPvalue_fc = new TLine(pvalue_fc_glob,pvalue_fc_glob,pvalue_fc_glob,pvalue_fc_glob);  
  TLatex *latZscore_fc = new TLatex(hdeltachisq_fc->GetXaxis()->GetBinLowEdge(5),hdeltachisq_fc->GetMaximum()*0.5, TString::Format("Z = %1.2f#sigma",zscore_fc_glob).Data() );
  latZscore_fc->Draw("same");
  latZscore_fc->SetTextSize(0.06);
  TLine *lZscore_fc = new TLine(zscore_fc_glob,zscore_fc_glob,zscore_fc_glob,zscore_fc_glob);
  TLegend *legPvalue_fc = new TLegend(0.5,0.6,0.8,0.8);
  legPvalue_fc->SetFillColor(0);
  legPvalue_fc->SetLineColor(0);
  legPvalue_fc->SetTextFont(132);
  legPvalue_fc->SetTextSize(0.06);
  sprintf(name,"#Delta #chi^{2}_{data} = %1.3f",deltachisq_fc_data); 
  legPvalue_fc->AddEntry(ldeltachisq_fc,name,"l");
  legPvalue_fc->Draw();  
  
  char outname[200];
  sprintf(outname,"%s.root",mytag.Data());
  TFile *outfile = new TFile(outname,"RECREATE");
  hdeltachisq_nh_ih->Write("hdeltachisq_Toys_nh_ih");
  ldeltachisq_nh_ih->Write("ldeltachisq_Data_nh_ih");
  lPvalue_nh_ih    ->Write("Pvalue_nh_ih");
  lZscore_nh_ih    ->Write("Zscore_nh_ih");
  hdeltachisq_fc   ->Write("hdeltachisq_Toys_fc");
  ldeltachisq_fc   ->Write("ldeltachisq_Data_fc");
  lPvalue_fc       ->Write("Pvalue_fc");
  lZscore_fc       ->Write("Zscore_fc");
  outfile->Close();

  // Print results
  
  cout << endl;
  cout << std::setprecision(5);
  cout << "Fit of: " << *trueparname0 << " Vs " << *trueparname1 << endl;
  cout << "- Toy experiments:" << endl;
  cout << "   Toy experiments are generated without marginalizing the Mass Hierarchy" << endl; 
  cout << "   True Mass Hierarchy: " << *trueMH << endl;
  cout << "   Total # of toy experiments = " << Ntoys << endl;
  cout << "   truepar0 = " << truepar0 << endl;
  cout << "   truepar1 = " << truepar1 << endl;
  cout << endl;
  cout << "- Data: " << endl;
  cout << "   Best-fit Mass Hierarchy: " << databfMH << endl;
  cout << "   DeltaChisq = Chisq(NH) - Chisq(IH) = " << deltachisq_nh_ih_data << endl;
  cout << "   DeltaChisq = Chisq(True) - Chisq(Best-Fit) = " << deltachisq_fc_data << endl;
  cout << endl;
  cout << "- Calculate right sided p-values:" << endl;
  cout << "   p-value = " << pvalue_nh_ih_right << endl;
  cout << "- Calculate left sided p-values:" << endl;
  cout << "   p-value = " << pvalue_nh_ih_left << endl;
  cout << "- The global p-value is calculated as ";
  if     (*trueMH == "normal")   cout << "right-side" << endl;
  else if(*trueMH == "inverted") cout << "left-side" << endl;
  cout << "- DeltaChisq = Chisq(NH) - Chisq(IH):" << endl;
  cout << "   p-value = " << pvalue_nh_ih_glob << endl;
  cout << "   z-score = " << zscore_nh_ih_glob << " sigma" << endl;
  cout << "- DeltaChisq = Chisq(True) - Chisq(Best-Fit):" << endl;
  cout << "   p-value = " << pvalue_fc_glob << endl;
  cout << "   z-score = " << zscore_fc_glob << " sigma" << endl;
  return;
}



double GetPvalue_rightside(vector<double> vec_toys,int ntoys,double data){
  double ntot_toys = ntoys;
  if(ntot_toys != vec_toys.size()){    
    cout << endl;
    cout << "Warning:" << endl;
    cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
    cout << "Input total number of toys = " << ntot_toys << endl; 
    ntot_toys = vec_toys.size();
    cout << "Set total number of toys to vector size: " << ntot_toys << endl;
    cout << endl;
  }
  if(ntot_toys<=0.){
    cerr << endl;
    cerr << "ntot_toys = " << ntot_toys << endl;
    cerr << endl;
    exit(1);
  }

  cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
  cout << "Input total number of toys = " << ntot_toys << endl; 
  ntot_toys = vec_toys.size();
  cout << "Set total number of toys to vector size: " << ntot_toys << endl;
  cout << endl;

  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it >= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}

double GetPvalue_leftside(vector<double> vec_toys,int ntoys,double data){
  double ntot_toys = ntoys;
  if(ntot_toys != vec_toys.size()){    
    cout << endl;
    cout << "Warning:" << endl;
    cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
    cout << "Input total number of toys = " << ntot_toys << endl; 
    ntot_toys = vec_toys.size();
    cout << "Set total number of toys to vector size: " << ntot_toys << endl;
    cout << endl;
  }
  if(ntot_toys<=0.){
    cerr << endl;
    cerr << "ntot_toys = " << ntot_toys << endl;
    cerr << endl;
    exit(1);
  }

  cout << "Total number of toys stored in vector = " << vec_toys.size() << endl; 
  cout << "Input total number of toys = " << ntot_toys << endl; 
  ntot_toys = vec_toys.size();
  cout << "Set total number of toys to vector size: " << ntot_toys << endl;
  cout << endl;

  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it <= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range 
  }
  double pvalue =  nfrac / (double)ntot_toys;

  return pvalue;
}

double GetPvalue_Discrete(TH1D *histo,double data){

  int bin_data = histo->FindBin(data);
  double pvalue = histo->Integral(bin_data,histo->GetNbinsX());
  pvalue /= histo->Integral();
  
  return pvalue;
}

