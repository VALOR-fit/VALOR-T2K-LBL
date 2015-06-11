//
// Sensitivity study w/ 1D or 2D plost
//


void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style);
TGraph * get_contour(TH2D * h);

void draw_sensitivity_plots_toys(
  const char *tag = "t23Vsdmsq_fakedataset_F_NuOsc",
  const bool PrintPlots = true,
  
  const bool NuBarMode = false,
 
  const bool doShiftMin = true,
  const double CLval = 4.61, //2.71, //4.61, // 90%CL
  const bool reactor = true, 
  const bool theta2theta = false,
  const char * parname1_def = "mixing13", //"mixing23",
  const char * parname2_def = "none", //"dmsq",
  const char * parname3_def = "none",
  
  const bool ChangeOrder_nh = false, // change order of params (2D plots)
  const bool ChangeOrder_ih = false, // change order of params (2D plots) 
  
  const bool DrawLeg = false,
  
  const char * LegendName_nh = "Fix NuOsc, Stat-Only",
  const char * LegendName_ih = "Fix NuOsc, Stat-Only", 
  
 
  const char * name_nh = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fakedataset_0p506_201411/fitgrid_NuOsc_0p506_fakedata_F_mixing23_dmsq_none_BinsTot_21_BinsTot_21_BinsTot_1_BinsPart_1_21_BinsPart_1_21_BinsPart_1_1",
  const char * name_ih = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/fakedataset_0p506_201411/fitgrid_NuOsc_0p506_fakedata_F_mixing23_dmsq_none_BinsTot_21_BinsTot_21_BinsTot_1_BinsPart_1_21_BinsPart_1_21_BinsPart_1_1",
 
  const char * output_file = "toys_sensitivity_plots.root"
  )
{
  using namespace std;
  using namespace std::string;
  using namespace valor;
  using namespace valor::conv;
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  string img_format(".eps");

  string filename_nh_s(Form("%s.root",name_nh));
  const char * const filename_nh = filename_nh_s.c_str();
  FitParamList * fparm_nh = 0;
  PhysParamList * iparm_nh = new PhysParamList;
  int nbins_param1_nh, nbins_param2_nh;
  int nbins_tot_param1_nh, nbins_tot_param2_nh;
  int nbins1_nh = 0;
  int nbins2_nh = 0;
  string *parname1_nh=new string();
  string *parname2_nh=new string();

  string filename_ih_s(Form("%s.root",name_ih));
  const char * const filename_ih = filename_ih_s.c_str();
  FitParamList * fparm_ih = 0;
  PhysParamList * iparm_ih = new PhysParamList;
  int nbins_param1_ih, nbins_param2_ih;
  int nbins_tot_param1_ih, nbins_tot_param2_ih;
  int nbins1_ih = 0;
  int nbins2_ih = 0;
  string *parname1_ih=new string();
  string *parname2_ih=new string();
  
  int max_bin_param1, min_bin_param1, max_bin_param2, min_bin_param2;
  double max_val_param1, min_val_param1, max_val_param2, min_val_param2;
  double step1, step2;
  
  PhysParamId_t param1,param2;
  
  param1 = util::PhysParamIdFromString(parname1_def);
  param2 = util::PhysParamIdFromString(parname2_def);
  
  if (NuBarMode) {
    param1 = util::GetBarredVersionOfPhysParam(param1);
    param2 = util::GetBarredVersionOfPhysParam(param2);
  }

  // 1D: dChi2 vs param1  (Normal Hierarchy)

  cout << endl;
  cout << "Read file: " << filename_nh << endl;
  cout << endl;
  
  TChain ch_nh("fitnt");

  if(ChangeOrder_nh){
    ch_nh.Add(filename_nh);
    ch_nh.SetBranchAddress("best_fit_params", &fparm_nh); 
    ch_nh.SetBranchAddress("nbins_param1", &nbins_param2_nh); 
    ch_nh.SetBranchAddress("nbins_param2", &nbins_param1_nh); 
    ch_nh.SetBranchAddress("nbins_tot_param1", &nbins_tot_param2_nh); 
    ch_nh.SetBranchAddress("nbins_tot_param2", &nbins_tot_param1_nh); 
    ch_nh.SetBranchAddress("parname1", &parname2_nh); // name of the grid parameter 2
    ch_nh.SetBranchAddress("parname2", &parname1_nh); // name of the grid parameter 1
    ch_nh.SetBranchAddress("max_val_param1",&max_val_param2);
    ch_nh.SetBranchAddress("min_val_param1",&min_val_param2);
    ch_nh.SetBranchAddress("max_val_param2",&max_val_param1);
    ch_nh.SetBranchAddress("min_val_param2",&min_val_param1);
    ch_nh.SetBranchAddress("step1",&step2);
    ch_nh.SetBranchAddress("step2",&step1);
  }
  else{
    ch_nh.Add(filename_nh);
    ch_nh.SetBranchAddress("best_fit_params", &fparm_nh); 
    ch_nh.SetBranchAddress("nbins_param1", &nbins_param1_nh); 
    ch_nh.SetBranchAddress("nbins_param2", &nbins_param2_nh); 
    ch_nh.SetBranchAddress("nbins_tot_param1", &nbins_tot_param1_nh); 
    ch_nh.SetBranchAddress("nbins_tot_param2", &nbins_tot_param2_nh); 
    ch_nh.SetBranchAddress("parname1", &parname1_nh); // name of the grid parameter 1
    ch_nh.SetBranchAddress("parname2", &parname2_nh); // name of the grid parameter 2
    ch_nh.SetBranchAddress("max_val_param1",&max_val_param1);
    ch_nh.SetBranchAddress("min_val_param1",&min_val_param1);
    ch_nh.SetBranchAddress("max_val_param2",&max_val_param2);
    ch_nh.SetBranchAddress("min_val_param2",&min_val_param2);
    ch_nh.SetBranchAddress("step1",&step1);
    ch_nh.SetBranchAddress("step2",&step2);
  }
  
  if(*parname1_nh==*parname2_nh && *parname1_nh!=""){
   cout << "Name parameter 1 == Name parameter 2 = " << *parname1_nh << endl; 
   exit(0);
  }

  const int nentries_nh = ch_nh.GetEntries();
  cout << "nentries (NH) = " << nentries_nh << endl;
  
  ch_nh.GetEntry(0);
  int nbins1_nh = TMath::Max(nbins_tot_param1_nh,nbins_param1_nh); // if using make_toygrid_toys_joint.C nbins_tot_param* is used instead of nbins_param*
  int nbins2_nh = TMath::Max(nbins_tot_param2_nh,nbins_param2_nh); // otherwise nbins_tot_param* = 1 < nbins_param*
  
  int nDim=0;
  if(nbins1_nh>1 || nbins2_nh>1) nDim = 1;
  if(nbins1_nh>1 && nbins2_nh>1) nDim = 2;
  
   if(*parname1_nh!=""){
    cout << endl;
    cout << "- Parameter 1: " << *parname1_nh << endl;
    cout << "  # bins = " << nbins1_nh << endl;
    cout << "  min = " << min_val_param1 << endl;
    cout << "  max = " << max_val_param1 << endl;
    cout << "  step = " << step1 << endl;
    cout << endl;
  } 
  if(*parname2_nh!=""){
  cout << "- Parameter 2: " << *parname2_nh << endl;
  cout << "  # bins = " << nbins2_nh << endl;
  cout << "  min = " << min_val_param2 << endl;
  cout << "  max = " << max_val_param2 << endl;
  cout << "  step = " << step2 << endl;
  cout << endl;
  }

  if(*parname1_nh=="" && *parname2_nh==""){
    cout << "# bins parameter 1 = " << nbins1_nh << endl; 
    cout << "# bins parameter 2 = " << nbins2_nh << endl; 
    cout << endl;
  }

  const int ntoys_nh = nentries_nh / (nbins1_nh * nbins2_nh);
  cout << "ntoys (NH) = " << ntoys_nh << endl;
  cout << endl;
  assert(ntoys_nh);
  
  //Set up bins and histograms
  double mixing13_min = 0;
  double mixing13_max = 0;
  double mixing13_step = 0;
  double mixing23_min = 0;
  double mixing23_max = 0;
  double mixing23_step = 0;
  double dmsq_min = 0;
  double dmsq_max = 0;
  double dmsq_step = 0;
  double deltacp_min = 0;
  double deltacp_max = 0;
  double deltacp_step = 0;

  if(!theta2theta){
    if(reactor){
      mixing13_min = 0.015;
      mixing13_max = 0.035;
      //mixing13_step = 0.001;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_nh-1);
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.1;
      //mixing13_step = 0.005;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_nh-1);
    }
  }
  else{
    if(reactor){
      mixing13_min = 0.06;
      mixing13_max = 0.14;
      //mixing13_step = 0.004;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_nh-1);
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.4;
      //mixing13_step = 0.02;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_nh-1);
    }
  }
  mixing23_min = 0.;
  mixing23_max = 1.;
  //mixing23_min = 0.3;
  //mixing23_max = 0.7;
  // //mixing23_step = 0.01;
  mixing23_step = (mixing23_max-mixing23_min)/(nbins1_nh-1);
  



  //dmsq_min = 2.0;//E-3;
  //dmsq_max = 3.0;//E-3;
  // //dmsq_step = 0.05;//E-3;

  dmsq_min = 3.0;//E-4;
  dmsq_max = 1.0;//E-4;
  //dmsq_step = 0.05;//E-2;
  dmsq_step = (dmsq_max-dmsq_min)/(nbins1_nh-1);
  
  deltacp_min = -conv::kpi;
  deltacp_max = conv::kpi;
  //deltacp_step = 2*conv::kpi/50;
  deltacp_step = (deltacp_max-deltacp_min)/(nbins1_nh-1);
  
  // Take steps and min-max values from the tree if they exist

  if(*parname1_nh==""){
    
    parname1_nh[0].append(parname1_def); // fill parname2 w/ default if not in the tree
    
    if(parname1_def=="deltacp"){
      min_val_param1 = deltacp_min;
      max_val_param1 = deltacp_max;
      step1 = deltacp_step;
    }
    else if(parname1_def=="mixing13"){
      min_val_param1 = mixing13_min;
      max_val_param1 = mixing13_max;
      step1 = mixing13_step;
    }
    else if(parname1_def=="mixing23"){
      min_val_param1 = mixing23_min;
      max_val_param1 = mixing23_max;
      step1 = mixing23_step;
    }
    else if(parname1_def=="dmsq"){
      min_val_param1 = dmsq_min;
      max_val_param1 = dmsq_max;
      step1 = dmsq_step;
    }
    else if(parname1_def=="none"){
      min_val_param1 = 0.;
      max_val_param1 = 0.;
      step1 = 0.;
    }
  }
  
  if(*parname2_nh==""){

    parname2_nh[0].append(parname2_def); // fill parname2 w/ default if not in the tree

    if(parname2_def=="deltacp"){
      min_val_param2 = deltacp_min;
      max_val_param2 = deltacp_max;
      step2 = deltacp_step;
    }
    else if(parname2_def=="mixing13"){
      min_val_param2 = mixing13_min;
      max_val_param2 = mixing13_max;
      step2 = mixing13_step;
    }
    else if(parname2_def=="mixing23"){
      min_val_param2 = mixing23_min;
      max_val_param2 = mixing23_max;
      step2 = mixing23_step;
    }
    else if(parname2_def=="dmsq"){
      min_val_param2 = dmsq_min;
      max_val_param2 = dmsq_max;
      step2 = dmsq_step;
    }
    else if(parname2_def=="none"){
    }
  }
   
  //
  char histoname[200];

  // 1D histograms

  Double_t** bestfit_points_nh_x = new Double_t*[ntoys_nh]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    bestfit_points_nh_x[itoy] = new Double_t[nbins1_nh];
  }
  Double_t** bestfit_points_nh_y = new Double_t*[ntoys_nh]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    bestfit_points_nh_y[itoy] = new Double_t[nbins1_nh];
  }

  sprintf(histoname,"h_chisq_%s_Ave_nh",parname1_nh[0].c_str());
  TH1D * h_chisq_par1_Ave_nh = new TH1D(histoname,histoname,nbins1_nh,min_val_param1-0.5*step1,max_val_param1+0.5*step1);

  TH1D ** h_chisq_par1_nh = new TH1D*[ntoys_nh];
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    sprintf(histoname,"h_%s_nh_Toy%d",parname1_nh[0].c_str(),itoy);
    h_chisq_par1_nh[itoy] = new TH1D(histoname,histoname,nbins1_nh,min_val_param1-0.5*step1,max_val_param1+0.5*step1);
  }  
  
  TH1D ** h_chisq_par1_nh_bin = new TH1D*[nbins1_nh];
  for(int ibin=0;ibin<nbins1_nh;ibin++){
    char histoname[200];
    sprintf(histoname,"h_chisq_%s_nh_bin%d",parname1_nh[0].c_str(),ibin);
    h_chisq_par1_nh_bin[ibin] = new TH1D(histoname,histoname,50,0,500);    
  }

  // 2D histograms

  sprintf(histoname,"h_chisq_%s_%s_Ave_nh",parname2_nh[0].c_str(),parname1_nh[0].c_str());
  TH2D * h_chisq_par2_par1_Ave_nh = new TH2D(histoname,histoname,nbins1_nh,min_val_param1-0.5*step1,max_val_param1+0.5*step1,nbins2_nh,min_val_param2-0.5*step2,max_val_param2+0.5*step2);

  TH2D ** h_chisq_par2_par1_nh= new TH2D*[ntoys_nh];
  for(int itoy=0; itoy<ntoys_nh; itoy++){
    char histoname[200];
    sprintf(histoname,"h_chisq_%s_%s_nh_Toy%d",parname2_nh[0].c_str(),parname1_nh[0].c_str(),itoy);
    h_chisq_par2_par1_nh[itoy] = new TH2D(histoname,histoname,nbins1_nh,min_val_param1-0.5*step1,max_val_param1+0.5*step1,nbins2_nh,min_val_param2-0.5*step2,max_val_param2+0.5*step2);
  }  

  TH1D *** h_chisq_par2_par1_nh_bin = new TH1D**[nbins1_nh];
  for(int ibin1=0;ibin1<nbins1_nh;ibin1++){
    h_chisq_par2_par1_nh_bin[ibin1] = new TH1D*[nbins2_nh];
    for(int ibin2=0;ibin2<nbins2_nh;ibin2++){
      char histoname[200];
      sprintf(histoname,"h_chisq_%s_%s_nh_bin%d_bin%d",parname2_nh[0].c_str(),parname1_nh[0].c_str(),ibin2,ibin1);
      h_chisq_par2_par1_nh_bin[ibin1][ibin2] = new TH1D(histoname,histoname,50,0,500);    
    }
  }

  // Loop over toys in the tree

  int iEvtToy = 0;
  for(int i = 0; i < nentries_nh; i++){
    iEvtToy = i / (nbins1_nh*nbins2_nh);
    if( !(i%(nbins1_nh*nbins2_nh)) ) cout << "# Toy " << iEvtToy << endl;

    ch_nh.GetEntry(i);
    int migradminstat_nh = fparm_nh->GetStatusCode("migrad_min");
    if (migradminstat_nh == 4) cout << "Fit failed!" << endl;
    //cout << "migradminstat_nh = " << migradminstat_nh << endl;

    double mix13_nh = fparm_nh->BestFitValue(conv::kPhypMixing13); // Parameters are fixed in the fit
    double mix23_nh = fparm_nh->BestFitValue(conv::kPhypMixing23);
    double dmsq_nh = fparm_nh->BestFitValue(conv::kPhypDm23sq);
    double deltacp_nh = fparm_nh->BestFitValue(conv::kPhypCPV);
    double chisq_nh = fparm_nh->GetBestFitStatistic();
    if(NuBarMode){
      mix13_nh = fparm_nh->BestFitValue(conv::kPhypMixing13Bar); // Parameters are fixed in the fit
      mix23_nh = fparm_nh->BestFitValue(conv::kPhypMixing23Bar);
      dmsq_nh = fparm_nh->BestFitValue(conv::kPhypDm23sqBar);
      deltacp_nh = fparm_nh->BestFitValue(conv::kPhypCPVBar);
      chisq_nh = fparm_nh->GetBestFitStatistic();
    }
    
 
    double bf1_nh=0;
    if(*parname1_nh=="deltacp") bf1_nh = deltacp_nh;
    else if(*parname1_nh=="mixing13") bf1_nh = mix13_nh;
    else if(*parname1_nh=="mixing23") bf1_nh = mix23_nh;
    else if(*parname1_nh=="dmsq") bf1_nh = dmsq_nh;
    double bf2_nh=0;
    if(*parname2_nh=="deltacp") bf2_nh = deltacp_nh;
    else if(*parname2_nh=="mixing13") bf2_nh = mix13_nh;
    else if(*parname2_nh=="mixing23") bf2_nh = mix23_nh;
    else if(*parname2_nh=="dmsq") bf2_nh = dmsq_nh;

    int bin_idx=0;
    double bincont=0.;
    if(nbins1_nh>1 && nbins2_nh<=1){ // 1D param1
      bin_idx = h_chisq_par1_nh[iEvtToy]->FindBin(bf1_nh);
      bincont = h_chisq_par1_nh[iEvtToy]->GetBinContent(bin_idx);
    }
    else if(nbins1_nh<=1 && nbins2_nh>1){ // 1D param2 
      bin_idx = h_chisq_par1_nh[iEvtToy]->FindBin(bf1_nh);
      bincont = h_chisq_par1_nh[iEvtToy]->GetBinContent(bin_idx);
    }
    else if(nbins1_nh>1 && nbins2_nh>1){ // 2D param2 Vs param1
      bin_idx = h_chisq_par2_par1_nh[iEvtToy]->FindBin(bf1_nh,bf2_nh);
      bincont = h_chisq_par2_par1_nh[iEvtToy]->GetBinContent(bin_idx);
    }
    if(bincont!=0){
      if(bincont>chisq_nh){
        if(chisq_nh!=0){
          if(nDim==1) h_chisq_par1_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
          else if(nDim==2) h_chisq_par2_par1_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
          //bestfit_points_nh_x[iEvtToy][bin_idx-1]  = bf1_nh;//sinsq2theta13;
          //bestfit_points_nh_y[iEvtToy][bin_idx-1]  = bf2_nh;
        }
      }
    } 
    else{ 
      if(chisq_nh!=0){
        if(nDim==1) h_chisq_par1_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
        else if(nDim==2) h_chisq_par2_par1_nh[iEvtToy]->SetBinContent(bin_idx,chisq_nh+bincont);
        //bestfit_points_nh_x[iEvtToy][bin_idx-1] = bf1_nh;//sinsq2theta13;
        //bestfit_points_nh_y[iEvtToy][bin_idx-1] = bf2_nh;
      }
    }
  } // end loop over toys
  
  
  // Get the dmsq scheme for each hierarchy
  NuOscDmsqScheme_t dmsqscheme_nh = kNuDmsqUndefined;
  if      (fparm_nh->Name(conv::kPhypDm23sq).find("32")!=string::npos) { dmsqscheme_nh = kNuDmsq12andAbs23;}
  else if (fparm_nh->Name(conv::kPhypDm23sq).find("31")!=string::npos) { dmsqscheme_nh = kNuDmsq12andAbs31;}
  else if (fparm_nh->Name(conv::kPhypDm23sq).find("/2")!=string::npos) { dmsqscheme_nh = kNuDmsqFogliLisi ;}
  else {cerr<<"Cannot determine NH mass splitting scheme";exit(-1);}

  
 
  // TGraph ** gbestfit_nh_Toy = new TGraph*[ntoys_nh]; 
  // for(int itoy=0;itoy<ntoys_nh;itoy++){
  //   char graphname[200];
  //   sprintf(graphname,"gbestfit_nh_Toy%d",itoy);
  //   gbestfit_nh_Toy[itoy] = new TGraph( nbins1_nh, bestfit_points_nh_x[itoy] ,bestfit_points_nh_y[itoy] ); // dcp Vs mix13  
  //   gbestfit_nh_Toy[itoy]->SetName(graphname);
  //   //cout << graphname << endl;
  // }  

  // Get the toy distribution for each param value
  if(nDim==1){
    cout << endl << endl;
    for(int ibin=0;ibin<nbins1_nh;ibin++){
      for(int itoy=0;itoy<ntoys_nh;itoy++){
        double bincont = h_chisq_par1_nh[itoy]->GetBinContent(ibin+1);
        h_chisq_par1_nh_bin[ibin]->Fill(bincont);
        cout << bincont << " ";
      }
    }
    cout << endl << endl;
  }
  if(nDim==2){
    for(int ibin1=0;ibin1<nbins1_nh;ibin1++){
      for(int ibin2=0;ibin2<nbins2_nh;ibin2++){
        for(int itoy=0;itoy<ntoys_nh;itoy++){
          int ibin = h_chisq_par2_par1_nh[itoy]->GetBin(ibin1+1,ibin2+1);
          double bincont = h_chisq_par2_par1_nh[itoy]->GetBinContent(ibin);
          h_chisq_par2_par1_nh_bin[ibin1][ibin2]->Fill(bincont);

          cout << bincont << endl;
        }
      }
    }
  }
  
  // Compute the average sensitivity    
  for(int itoy=0;itoy<ntoys_nh;itoy++){
    h_chisq_par1_Ave_nh->Add(h_chisq_par1_nh[itoy]);
    h_chisq_par2_par1_Ave_nh->Add(h_chisq_par2_par1_nh[itoy]);
  }
  h_chisq_par1_Ave_nh->Scale(1./ntoys_nh);
  h_chisq_par2_par1_Ave_nh->Scale(1./ntoys_nh);
  
  
  // Set plot title
  {
  NuOscDmsqScheme_t orig_scheme = Opt::Instance()->gOscDmsqScheme; // save the scheme so we can put it back later
  Opt::Instance()->gOscDmsqScheme = dmsqscheme_nh;
  h_chisq_par1_Ave_nh->GetXaxis()->SetTitle(util::AsStringPaperRooTex(param1).c_str());
  Opt::Instance()->gOscDmsqScheme = orig_scheme;
  }
  // Find the minimum chi2 (global scan)
  
  cout << endl;
  cout << "Find the minimum chi2 (global scan)" << endl;
  cout << endl;

  double *min_chisq_nh = new double[ntoys_nh];
  double min_chisq_Ave_nh=999.;
    
  if(nDim==1){ 
    for(int itoy=0;itoy<ntoys_nh;itoy++){
      min_chisq_nh[itoy] = 999.;
      //cout << "# Toy " << itoy << endl;    
      for(int i = 1; i <= nbins1_nh; i++){
        double bin = h_chisq_par1_nh[itoy]->GetBinContent(i);
        if (bin < min_chisq_nh[itoy]) min_chisq_nh[itoy] = bin;
        //cout << bin << endl;
      }
      //cout << "For param1, min chisq  = " << min_chisq_nh[itoy] << endl;
      //cout << endl;
    }
    for(int i = 1; i <= nbins1_nh; i++){
      double bin = h_chisq_par1_Ave_nh->GetBinContent(i);
      if (bin < min_chisq_Ave_nh) min_chisq_Ave_nh = bin;
    }
    cout << "Global scanned chi2 minimum = " << min_chisq_Ave_nh << endl;     
    
  }
  else if(nDim==2){
    for(int itoy=0;itoy<ntoys_nh;itoy++){
      min_chisq_nh[itoy] = 999.;
      //cout << "# Toy " << itoy << endl;    
      for(int i1 = 1; i1 <= nbins1_nh; i1++){
        for(int i2 = 1; i2 <= nbins2_nh; i2++){
          double bin = h_chisq_par2_par1_nh[itoy]->GetBinContent(i1,i2);
          if (bin < min_chisq_nh[itoy]) min_chisq_nh[itoy] = bin;
          //cout << bin << endl;
        }
      }
      //cout << "For param1 Vs param2, min chisq  = " << min_chisq_nh[itoy] << endl;
      //cout << endl;
    }

    for(int i1 = 1; i1 <= nbins1_nh; i1++){
      for(int i2 = 1; i2 <= nbins2_nh; i2++){
        double bin = h_chisq_par2_par1_Ave_nh->GetBinContent(i1,i2);
        if (bin < min_chisq_Ave_nh) min_chisq_Ave_nh = bin;
        cout << bin << " ";
      }
      cout << endl;
    }
    cout << "Global scanned chi2 minimum = " << min_chisq_Ave_nh << endl; 
  } // nDim==2
  

  // Take RMS of dcp in each true bin
  // Chi2 distribution is gaussian (many bins)

  TCanvas **c_chisq_par1_nh_bin = new TCanvas*[nbins1_nh];
  double *RMS_nh = new double[nbins1_nh];
  double *RelErr_nh = new double[nbins1_nh];      

  if(nDim==1){
    for(int ibin=0;ibin<nbins1_nh;ibin++){
      //~ char canvname[200];
      //~ sprintf(canvname,"c_chisq_%s_nh_bin%d",parname1_nh[0].c_str(),ibin);
      RMS_nh[ibin] = h_chisq_par1_nh_bin[ibin]->GetRMS();
      RelErr_nh[ibin] = h_chisq_par1_nh_bin[ibin]->GetRMS() / h_chisq_par1_nh_bin[ibin]->GetMean();
    }
  }

  // Sensi with error bars

  TGraphErrors *g_chisq_par1_nh_Err = new TGraphErrors();
  TGraph *g_chisq_par1_nh = new TGraph();
  if(nDim==1){
    for(int ipt=0;ipt<nbins1_nh;ipt++){
      double contx = h_chisq_par1_Ave_nh->GetBinCenter(ipt+1);
      double conty = h_chisq_par1_Ave_nh->GetBinContent(ipt+1);
      if(doShiftMin) conty -= min_chisq_Ave_nh; // shift the chi2 wrt the minimum
      double errx = h_chisq_par1_Ave_nh->GetBinWidth(ipt+1)/2;
      double erry = RMS_nh[ipt];
      
      //cout << contx << " +- " << errx << " - " << conty << " +- " << erry << endl;
      
      if(erry<=0.) erry = 0.000001;
      g_chisq_par1_nh->SetPoint(ipt,contx,conty);
      g_chisq_par1_nh_Err->SetPoint(ipt,contx,conty);
      g_chisq_par1_nh_Err->SetPointError(ipt,errx,erry);	
      //cout << contx << " +- " << errx << endl;
    }
  }
  
  // Shift the chi2 wrt the minimum
  if(doShiftMin){

    if(nDim==1){
      for(int itoy=0;itoy<ntoys_nh;itoy++){
        for(int i = 1; i <= nbins1_nh; i++){
          double bin = h_chisq_par1_nh[itoy]->GetBinContent(i);
          h_chisq_par1_nh[itoy]->SetBinContent(i,bin-min_chisq_nh[itoy]); // scanned global minimum
          //h_chisq_par1_nh[itoy]->SetBinContent(i,bin-chisq_bf_nh); // fiited global minimum
        }
      }
      for(int i = 1; i <= nbins1_nh; i++){
        double bin = h_chisq_par1_Ave_nh->GetBinContent(i);
        h_chisq_par1_Ave_nh->SetBinContent(i,bin-min_chisq_Ave_nh); // scanned global minimum
      }
    }
   
    else if(nDim==2){
      for(int itoy=0;itoy<ntoys_nh;itoy++){
        for(int i1 = 1; i1 <= nbins1_nh; i1++){
          for(int i2 = 1; i2 <= nbins2_nh; i2++){
            double bin = h_chisq_par2_par1_nh[itoy]->GetBinContent(i1,i2);
            h_chisq_par2_par1_nh[itoy]->SetBinContent(i1,i2,bin-min_chisq_nh[itoy]); // scanned global minimum
            //h_chisq_par1_nh[itoy]->SetBinContent(i,bin-chisq_bf_nh); // fiited global minimum
          }
        }
      }
      for(int i1 = 1; i1 <= nbins1_nh; i1++){
        for(int i2 = 1; i2 <= nbins2_nh; i2++){
          double bin = h_chisq_par2_par1_Ave_nh->GetBinContent(i1,i2);
          h_chisq_par2_par1_Ave_nh->SetBinContent(i1,i2,bin-min_chisq_Ave_nh); // scanned global minimum
        }
      }
    }
  }

  // 1D: dChi2 vs param1  (Inverted Hierarchy)

  cout << endl;
  cout << "Read file: " << filename_ih << endl;
  cout << endl;
  
  TChain ch_ih("fitnt");

  if(ChangeOrder_ih){
    ch_ih.Add(filename_ih);
    ch_ih.SetBranchAddress("best_fit_params", &fparm_ih); 
    ch_ih.SetBranchAddress("nbins_param1", &nbins_param2_ih); 
    ch_ih.SetBranchAddress("nbins_param2", &nbins_param1_ih); 
    ch_ih.SetBranchAddress("nbins_tot_param1", &nbins_tot_param2_ih); 
    ch_ih.SetBranchAddress("nbins_tot_param2", &nbins_tot_param1_ih); 
    ch_ih.SetBranchAddress("best_fit_params", &fparm_ih); 
    ch_ih.SetBranchAddress("parname1", &parname2_ih); // name of the grid parameter 2
    ch_ih.SetBranchAddress("parname2", &parname1_ih); // name of the grid parameter 1
    ch_ih.SetBranchAddress("max_val_param1",&max_val_param2);
    ch_ih.SetBranchAddress("min_val_param1",&min_val_param2);
    ch_ih.SetBranchAddress("max_val_param2",&max_val_param1);
    ch_ih.SetBranchAddress("min_val_param2",&min_val_param1);
    ch_ih.SetBranchAddress("step1",&step2);
    ch_ih.SetBranchAddress("step2",&step1);
  }
  else{
    ch_ih.Add(filename_ih);
    ch_ih.SetBranchAddress("best_fit_params", &fparm_ih); 
    ch_ih.SetBranchAddress("nbins_param1", &nbins_param1_ih); 
    ch_ih.SetBranchAddress("nbins_param2", &nbins_param2_ih); 
    ch_ih.SetBranchAddress("nbins_tot_param1", &nbins_tot_param1_ih); 
    ch_ih.SetBranchAddress("nbins_tot_param2", &nbins_tot_param2_ih); 
    ch_ih.SetBranchAddress("best_fit_params", &fparm_ih); 
    ch_ih.SetBranchAddress("parname1", &parname1_ih); // name of the grid parameter 1
    ch_ih.SetBranchAddress("parname2", &parname2_ih); // name of the grid parameter 2
    ch_ih.SetBranchAddress("max_val_param1",&max_val_param1);
    ch_ih.SetBranchAddress("min_val_param1",&min_val_param1);
    ch_ih.SetBranchAddress("max_val_param2",&max_val_param2);
    ch_ih.SetBranchAddress("min_val_param2",&min_val_param2);
    ch_ih.SetBranchAddress("step1",&step1);
    ch_ih.SetBranchAddress("step2",&step2);
  }
  
  if(*parname1_ih==*parname2_ih && *parname1_ih!=""){
   cout << "Name parameter 1 == Name parameter 2 = " << *parname1_ih << endl; 
   exit(0);
  }

  const int nentries_ih = ch_ih.GetEntries();
  cout << "nentries (IH) = " << nentries_ih << endl;
  
  ch_ih.GetEntry(0);
  int nbins1_ih = TMath::Max(nbins_tot_param1_ih,nbins_param1_ih); // if using make_toygrid_toys_joint.C nbins_tot_param* is used instead of nbins_param*
  int nbins2_ih = TMath::Max(nbins_tot_param2_ih,nbins_param2_ih); // otherwise nbins_tot_param* = 1 < nbins_param*
  
  int nDim=0;
  if(nbins1_ih>1 || nbins2_ih>1) nDim = 1;
  if(nbins1_ih>1 && nbins2_ih>1) nDim = 2;
 
  if(*parname1_ih!=""){
    cout << endl;
    cout << "- Parameter 1: " << *parname1_ih << endl;
    cout << "  # bins = " << nbins1_ih << endl;
    cout << "  min = " << min_val_param1 << endl;
    cout << "  max = " << max_val_param1 << endl;
    cout << "  step = " << step1 << endl;
    cout << endl;
  } 
  if(*parname2_ih!=""){
  cout << "- Parameter 2: " << *parname2_ih << endl;
  cout << "  # bins = " << nbins2_ih << endl;
  cout << "  min = " << min_val_param2 << endl;
  cout << "  max = " << max_val_param2 << endl;
  cout << "  step = " << step2 << endl;
  cout << endl;
  }

  if(*parname1_ih=="" && *parname2_ih==""){
    cout << "# bins parameter 1 = " << nbins1_ih << endl; 
    cout << "# bins parameter 2 = " << nbins2_ih << endl; 
    cout << endl;
  }

  const int ntoys_ih = nentries_ih / (nbins1_ih * nbins2_ih);
  cout << "ntoys (IH) = " << ntoys_ih << endl;
  cout << endl;
  assert(ntoys_ih);
  
  //Set up bins and histograms
  double mixing13_min = 0;
  double mixing13_max = 0;
  double mixing13_step = 0;
  double mixing23_min = 0;
  double mixing23_max = 0;
  double mixing23_step = 0;
  double dmsq_min = 0;
  double dmsq_max = 0;
  double dmsq_step = 0;
  double deltacp_min = 0;
  double deltacp_max = 0;
  double deltacp_step = 0;

  if(!theta2theta){
    if(reactor){
      mixing13_min = 0.015;
      mixing13_max = 0.035;
      //mixing13_step = 0.001;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_ih-1);
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.1;
      //mixing13_step = 0.005;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_ih-1);
    }
  }
  else{
    if(reactor){
      mixing13_min = 0.06;
      mixing13_max = 0.14;
      //mixing13_step = 0.004;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_ih-1);
    }
    else{
      mixing13_min = 0.0;
      mixing13_max = 0.4;
      //mixing13_step = 0.02;
      mixing13_step = (mixing13_max-mixing13_min)/(nbins1_ih-1);
    }
  }
  mixing23_min = 0.3;
  mixing23_max = 0.7;
  //mixing23_step = 0.01;
  mixing23_step = (mixing23_max-mixing23_min)/(nbins1_ih-1);
  
  dmsq_min = 2.0;//E-3;
  dmsq_max = 3.0;//E-3;
  //dmsq_step = 0.05;//E-3;
  dmsq_step = (dmsq_max-dmsq_min)/(nbins1_ih-1);

  deltacp_min = -conv::kpi;
  deltacp_max = conv::kpi;
  //deltacp_step = 2*conv::kpi/50;
  deltacp_step = (deltacp_max-deltacp_min)/(nbins1_ih-1);
  

  // Take steps and min-max values from the tree if they exist

  if(*parname1_ih==""){
    
    parname1_ih[0].append(parname1_def); // fill parname2 w/ default if not in the tree
    
    if(parname1_def=="deltacp"){
      min_val_param1 = deltacp_min;
      max_val_param1 = deltacp_max;
      step1 = deltacp_step;
    }
    else if(parname1_def=="mixing13"){
      min_val_param1 = mixing13_min;
      max_val_param1 = mixing13_max;
      step1 = mixing13_step;
    }
    else if(parname1_def=="mixing23"){
      min_val_param1 = mixing23_min;
      max_val_param1 = mixing23_max;
      step1 = mixing23_step;
    }
    else if(parname1_def=="dmsq"){
      min_val_param1 = dmsq_min;
      max_val_param1 = dmsq_max;
      step1 = dmsq_step;
    }
    else if(parname1_def=="none"){
      min_val_param1 = 0.;
      max_val_param1 = 0.;
      step1 = 0.;
    }
  }
  
  if(*parname2_ih==""){

    parname2_ih[0].append(parname2_def); // fill parname2 w/ default if not in the tree

    if(parname2_def=="deltacp"){
      min_val_param2 = deltacp_min;
      max_val_param2 = deltacp_max;
      step2 = deltacp_step;
    }
    else if(parname2_def=="mixing13"){
      min_val_param2 = mixing13_min;
      max_val_param2 = mixing13_max;
      step2 = mixing13_step;
    }
    else if(parname2_def=="mixing23"){
      min_val_param2 = mixing23_min;
      max_val_param2 = mixing23_max;
      step2 = mixing23_step;
    }
    else if(parname2_def=="dmsq"){
      min_val_param2 = dmsq_min;
      max_val_param2 = dmsq_max;
      step2 = dmsq_step;
    }
    else if(parname2_def=="none"){
    }
  }
    
  //
  char histoname[200];

  // 1D histograms

  Double_t** bestfit_points_ih_x = new Double_t*[ntoys_ih]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_ih; itoy++){
    bestfit_points_ih_x[itoy] = new Double_t[nbins1_ih];
  }
  Double_t** bestfit_points_ih_y = new Double_t*[ntoys_ih]; // bestfit points for each toy
  for(int itoy=0; itoy<ntoys_ih; itoy++){
    bestfit_points_ih_y[itoy] = new Double_t[nbins1_ih];
  }

  sprintf(histoname,"h_chisq_%s_Ave_ih",parname1_ih[0].c_str());
  TH1D * h_chisq_par1_Ave_ih = new TH1D(histoname,histoname,nbins1_ih,min_val_param1-0.5*step1,max_val_param1+0.5*step1);

  TH1D ** h_chisq_par1_ih = new TH1D*[ntoys_ih];
  for(int itoy=0; itoy<ntoys_ih; itoy++){
    sprintf(histoname,"h_%s_ih_Toy%d",parname1_ih[0].c_str(),itoy);
    h_chisq_par1_ih[itoy] = new TH1D(histoname,histoname,nbins1_ih,min_val_param1-0.5*step1,max_val_param1+0.5*step1);
  }  
  
  TH1D ** h_chisq_par1_ih_bin = new TH1D*[nbins1_ih];
  for(int ibin=0;ibin<nbins1_ih;ibin++){
    char histoname[200];
    sprintf(histoname,"h_chisq_%s_ih_bin%d",parname1_ih[0].c_str(),ibin);
    h_chisq_par1_ih_bin[ibin] = new TH1D(histoname,histoname,50,0,500);    
  }

  // 2D histograms

  sprintf(histoname,"h_chisq_%s_%s_Ave_ih",parname2_ih[0].c_str(),parname1_ih[0].c_str());
  TH2D * h_chisq_par2_par1_Ave_ih = new TH2D(histoname,histoname,nbins1_ih,min_val_param1-0.5*step1,max_val_param1+0.5*step1,nbins2_ih,min_val_param2-0.5*step2,max_val_param2+0.5*step2);

  TH2D ** h_chisq_par2_par1_ih= new TH2D*[ntoys_ih];
  for(int itoy=0; itoy<ntoys_ih; itoy++){
    char histoname[200];
    sprintf(histoname,"h_chisq_%s_%s_ih_Toy%d",parname2_ih[0].c_str(),parname1_ih[0].c_str(),itoy);
    h_chisq_par2_par1_ih[itoy] = new TH2D(histoname,histoname,nbins1_ih,min_val_param1-0.5*step1,max_val_param1+0.5*step1,nbins2_ih,min_val_param2-0.5*step2,max_val_param2+0.5*step2);
  }  

  TH1D *** h_chisq_par2_par1_ih_bin = new TH1D**[nbins1_ih];
  for(int ibin1=0;ibin1<nbins1_ih;ibin1++){
    h_chisq_par2_par1_ih_bin[ibin1] = new TH1D*[nbins2_ih];
    for(int ibin2=0;ibin2<nbins2_ih;ibin2++){
      char histoname[200];
      sprintf(histoname,"h_chisq_%s_%s_ih_bin%d_bin%d",parname2_ih[0].c_str(),parname1_ih[0].c_str(),ibin2,ibin1);
      h_chisq_par2_par1_ih_bin[ibin1][ibin2] = new TH1D(histoname,histoname,50,0,500);    
    }
  }

  // Loop over toys in the tree

  int iEvtToy = 0;
  for(int i = 0; i < nentries_ih; i++){
    iEvtToy = i / (nbins1_ih*nbins2_ih);
    if( !(i%(nbins1_ih*nbins2_ih)) ) cout << "# Toy " << iEvtToy << endl;
        
    ch_ih.GetEntry(i);
    int migradminstat_ih = fparm_ih->GetStatusCode("migrad_min");
    if (migradminstat_ih == 4) cout << "Fit failed!" << endl;
    //cout << "migradminstat_ih = " << migradminstat_ih << endl;

    double mix13_ih = fparm_ih->BestFitValue(conv::kPhypMixing13); // Parameters are fixed in the fit
    double mix23_ih = fparm_ih->BestFitValue(conv::kPhypMixing23);
    double dmsq_ih = fparm_ih->BestFitValue(conv::kPhypDm23sq);
    double deltacp_ih = fparm_ih->BestFitValue(conv::kPhypCPV);
    double chisq_ih = fparm_ih->GetBestFitStatistic();
    if(NuBarMode){
      mix13_ih = fparm_ih->BestFitValue(conv::kPhypMixing13Bar); // Parameters are fixed in the fit
      mix23_ih = fparm_ih->BestFitValue(conv::kPhypMixing23Bar);
      dmsq_ih = fparm_ih->BestFitValue(conv::kPhypDm23sqBar);
      deltacp_ih = fparm_ih->BestFitValue(conv::kPhypCPVBar);
      chisq_ih = fparm_ih->GetBestFitStatistic();
    }

    double bf1_ih=0;
    if(*parname1_ih=="deltacp") bf1_ih = deltacp_ih;
    else if(*parname1_ih=="mixing13") bf1_ih = mix13_ih;
    else if(*parname1_ih=="mixing23") bf1_ih = mix23_ih;
    else if(*parname1_ih=="dmsq") bf1_ih = dmsq_ih;
    double bf2_ih=0;
    if(*parname2_ih=="deltacp") bf2_ih = deltacp_ih;
    else if(*parname2_ih=="mixing13") bf2_ih = mix13_ih;
    else if(*parname2_ih=="mixing23") bf2_ih = mix23_ih;
    else if(*parname2_ih=="dmsq") bf2_ih = dmsq_ih;

    int bin_idx=0;
    double bincont=0.;

    if(nbins1_ih>1 && nbins2_ih<=1){ // 1D param1
      bin_idx = h_chisq_par1_ih[iEvtToy]->FindBin(bf1_ih);
      bincont = h_chisq_par1_ih[iEvtToy]->GetBinContent(bin_idx);
    }
    else if(nbins1_ih<=1 && nbins2_ih>1){ // 1D param2 
      bin_idx = h_chisq_par1_ih[iEvtToy]->FindBin(bf1_ih);
      bincont = h_chisq_par1_ih[iEvtToy]->GetBinContent(bin_idx);
    }
    else if(nbins1_ih>1 && nbins2_ih>1){ // 2D param2 Vs param1
      bin_idx = h_chisq_par2_par1_ih[iEvtToy]->FindBin(bf1_ih,bf2_ih);
      bincont = h_chisq_par2_par1_ih[iEvtToy]->GetBinContent(bin_idx);
    }
  
    if(bincont!=0){
      if(bincont>chisq_ih){
        if(chisq_ih!=0){
          if(nDim==1) h_chisq_par1_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
          else if(nDim==2) h_chisq_par2_par1_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
          //bestfit_points_ih_x[iEvtToy][bin_idx-1]  = bf1_ih;//sinsq2theta13;
          //bestfit_points_ih_y[iEvtToy][bin_idx-1]  = bf2_ih;
        }
      }
    } 
    else{ 
      if(chisq_ih!=0){
        if(nDim==1) h_chisq_par1_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
        else if(nDim==2) h_chisq_par2_par1_ih[iEvtToy]->SetBinContent(bin_idx,chisq_ih+bincont);
        //bestfit_points_ih_x[iEvtToy][bin_idx-1] = bf1_ih;//sinsq2theta13;
        //bestfit_points_ih_y[iEvtToy][bin_idx-1] = bf2_ih;
      }
    }
  } // end loop over toys
  
  // Get IH DMsq scheme
  NuOscDmsqScheme_t dmsqscheme_ih = kNuDmsqUndefined;
  if      (fparm_ih->Name(conv::kPhypDm23sq).find("32")!=string::npos) { dmsqscheme_ih = kNuDmsq12andAbs23;}
  else if (fparm_ih->Name(conv::kPhypDm23sq).find("31")!=string::npos) { dmsqscheme_ih = kNuDmsq12andAbs31;}
  else if (fparm_ih->Name(conv::kPhypDm23sq).find("/2")!=string::npos) { dmsqscheme_ih = kNuDmsqFogliLisi ;}
  else {cerr<<"Cannot determine IH mass splitting scheme";exit(-1);}
  
 
  // TGraph ** gbestfit_ih_Toy = new TGraph*[ntoys_ih]; 
  // for(int itoy=0;itoy<ntoys_ih;itoy++){
  //   char graphname[200];
  //   sprintf(graphname,"gbestfit_ih_Toy%d",itoy);
  //   gbestfit_ih_Toy[itoy] = new TGraph( nbins1_ih, bestfit_points_ih_x[itoy] ,bestfit_points_ih_y[itoy] ); // dcp Vs mix13  
  //   gbestfit_ih_Toy[itoy]->SetName(graphname);
  //   //cout << graphname << endl;
  // }  

  // Get the toy distribution for each param value
  if(nDim==1){
    cout << endl << endl;
    for(int ibin=0;ibin<nbins1_ih;ibin++){
      for(int itoy=0;itoy<ntoys_ih;itoy++){
        double bincont = h_chisq_par1_ih[itoy]->GetBinContent(ibin+1);
        h_chisq_par1_ih_bin[ibin]->Fill(bincont);

        cout << bincont << " ";
      }
    }
    cout << endl << endl;
  }
  if(nDim==2){
    for(int ibin1=0;ibin1<nbins1_ih;ibin1++){
      for(int ibin2=0;ibin2<nbins2_ih;ibin2++){
        for(int itoy=0;itoy<ntoys_ih;itoy++){
          int ibin = h_chisq_par2_par1_ih[itoy]->GetBin(ibin1+1,ibin2+1);
          double bincont = h_chisq_par2_par1_ih[itoy]->GetBinContent(ibin);
          h_chisq_par2_par1_ih_bin[ibin1][ibin2]->Fill(bincont);
        }
      }
    }
  }
  
  // Compute the average sensitivity    
  for(int itoy=0;itoy<ntoys_ih;itoy++){
    h_chisq_par1_Ave_ih->Add(h_chisq_par1_ih[itoy]);
    h_chisq_par2_par1_Ave_ih->Add(h_chisq_par2_par1_ih[itoy]);
  }
  h_chisq_par1_Ave_ih->Scale(1./ntoys_ih);
  h_chisq_par2_par1_Ave_ih->Scale(1./ntoys_ih);
  
  h_chisq_par1_Ave_ih->GetXaxis()->SetTitle(util::AsStringPaperRooTex(param1).c_str());
  
  // Find the minimum chi2 (global scan)
  
  cout << endl;
  cout << "Find the minimum chi2 (global scan)" << endl;
  cout << endl;

  double *min_chisq_ih = new double[ntoys_ih];
  double min_chisq_Ave_ih=999.;
    
  if(nDim==1){ 
    for(int itoy=0;itoy<ntoys_ih;itoy++){
      min_chisq_ih[itoy] = 999.;
      //cout << "# Toy " << itoy << endl;    
      for(int i = 1; i <= nbins1_ih; i++){
        double bin = h_chisq_par1_ih[itoy]->GetBinContent(i);
        if (bin < min_chisq_ih[itoy]) min_chisq_ih[itoy] = bin;
        //cout << bin << endl;
      }
      //cout << "For param1, min chisq  = " << min_chisq_ih[itoy] << endl;
      //cout << endl;
    }
    for(int i = 1; i <= nbins1_ih; i++){
      double bin = h_chisq_par1_Ave_ih->GetBinContent(i);
      if (bin < min_chisq_Ave_ih) min_chisq_Ave_ih = bin;
    }
    cout << "Global scanned chi2 minimum = " << min_chisq_Ave_ih << endl;     
    
    // Set plot title
    {
    NuOscDmsqScheme_t orig_scheme = Opt::Instance()->gOscDmsqScheme; // save the scheme so we can put it back later
    Opt::Instance()->gOscDmsqScheme = dmsqscheme_ih;
    h_chisq_par1_Ave_ih->GetXaxis()->SetTitle(util::AsStringPaperRooTex(param1).c_str());
    Opt::Instance()->gOscDmsqScheme = orig_scheme;
    }
    
  } // nDim==1
  else if(nDim==2){
    for(int itoy=0;itoy<ntoys_ih;itoy++){
      min_chisq_ih[itoy] = 999.;
      //cout << "# Toy " << itoy << endl;    
      for(int i1 = 1; i1 <= nbins1_ih; i1++){
        for(int i2 = 1; i2 <= nbins2_ih; i2++){
          double bin = h_chisq_par2_par1_ih[itoy]->GetBinContent(i1,i2);
          if (bin < min_chisq_ih[itoy]) min_chisq_ih[itoy] = bin;
          //cout << bin << endl;
        }
      }
      //cout << "For param1 Vs param2, min chisq  = " << min_chisq_ih[itoy] << endl;
      //cout << endl;
    }

    for(int i1 = 1; i1 <= nbins1_ih; i1++){
      for(int i2 = 1; i2 <= nbins2_ih; i2++){
        double bin = h_chisq_par2_par1_Ave_ih->GetBinContent(i1,i2);
        if (bin < min_chisq_Ave_ih) min_chisq_Ave_ih = bin;
        cout << bin << " ";
      }
      cout << endl;
    }
    cout << "Global scanned chi2 minimum = " << min_chisq_Ave_ih << endl; 
  } // nDim==2
  
  // Take RMS of dcp in each true bin
  // Chi2 distribution is gaussian (many bins)

  TCanvas **c_chisq_par1_ih_bin = new TCanvas*[nbins1_ih];
  double *RMS_ih = new double[nbins1_ih];
  double *RelErr_ih = new double[nbins1_ih];      

  if(nDim==1){
    for(int ibin=0;ibin<nbins1_ih;ibin++){
      //~ char canvname[200];
      //~ sprintf(canvname,"c_chisq_%s_ih_bin%d",parname1_ih[0].c_str(),ibin);
      RMS_ih[ibin] = h_chisq_par1_ih_bin[ibin]->GetRMS();
      RelErr_ih[ibin] = h_chisq_par1_ih_bin[ibin]->GetRMS() / h_chisq_par1_ih_bin[ibin]->GetMean();
    }
  }

  // Sensi with error bars

  TGraphErrors *g_chisq_par1_ih_Err = new TGraphErrors();
  TGraph *g_chisq_par1_ih = new TGraph();
  if(nDim==1){
    for(int ipt=0;ipt<nbins1_ih;ipt++){
      double contx = h_chisq_par1_Ave_ih->GetBinCenter(ipt+1);
      double conty = h_chisq_par1_Ave_ih->GetBinContent(ipt+1);
      if(doShiftMin) conty -= min_chisq_Ave_ih; // shift the chi2 wrt the minimum
      double errx = h_chisq_par1_Ave_ih->GetBinWidth(ipt+1)/2;
      double erry = RMS_ih[ipt];
      
      //cout << contx << " +- " << errx << " - " << conty << " +- " << erry << endl;
      
      if(erry<=0.) erry = 0.000001;
      g_chisq_par1_ih->SetPoint(ipt,contx,conty);
      g_chisq_par1_ih_Err->SetPoint(ipt,contx,conty);
      g_chisq_par1_ih_Err->SetPointError(ipt,errx,erry);	
      //cout << contx << " +- " << errx << endl;
    }
  }
  
  //Shift the chi2 wrt the minimum
  if(doShiftMin){
    if(nDim==1){
      for(int itoy=0;itoy<ntoys_ih;itoy++){
        for(int i = 1; i <= nbins1_ih; i++){
          double bin = h_chisq_par1_ih[itoy]->GetBinContent(i);
          h_chisq_par1_ih[itoy]->SetBinContent(i,bin-min_chisq_ih[itoy]); // scanned global minimum
          //h_chisq_par1_ih[itoy]->SetBinContent(i,bin-chisq_bf_ih); // fitted global minimum
        }
      }
      for(int i = 1; i <= nbins1_ih; i++){
        double bin = h_chisq_par1_Ave_ih->GetBinContent(i);
        h_chisq_par1_Ave_ih->SetBinContent(i,bin-min_chisq_Ave_ih); // scanned global minimum
      }
    }
    else if(nDim==2){
      for(int itoy=0;itoy<ntoys_ih;itoy++){
        for(int i1 = 1; i1 <= nbins1_ih; i1++){
          for(int i2 = 1; i2 <= nbins2_ih; i2++){
            double bin = h_chisq_par2_par1_ih[itoy]->GetBinContent(i1,i2);
            h_chisq_par2_par1_ih[itoy]->SetBinContent(i1,i2,bin-min_chisq_ih[itoy]); // scanned global minimum
            //h_chisq_par1_ih[itoy]->SetBinContent(i,bin-chisq_bf_ih); // fitted global minimum
          }
        }
      }
      for(int i1 = 1; i1 <= nbins1_ih; i1++){
        for(int i2 = 1; i2 <= nbins2_ih; i2++){
          double bin = h_chisq_par2_par1_Ave_ih->GetBinContent(i1,i2);
          h_chisq_par2_par1_Ave_ih->SetBinContent(i1,i2,bin-min_chisq_Ave_ih); // scanned global minimum
        }
      }
    }
  }
  

  // Shift between NH and IH
  TGraph *g_chisq_par1_ih_shift = new TGraph();
  TGraph *g_chisq_par1_nh_shift = new TGraph();
  if(doShiftMin && nDim==1){
    double shift_nh=0.;
    double shift_ih=0.;
    if(min_chisq_Ave_ih<min_chisq_Ave_nh) shift_nh = min_chisq_Ave_nh - min_chisq_Ave_ih;
    else                                  shift_ih = min_chisq_Ave_ih - min_chisq_Ave_nh;
    for(int ipt=0;ipt<nbins1_nh;ipt++){      
      double valx = g_chisq_par1_nh->GetX()[ipt];
      double valy = g_chisq_par1_nh->GetY()[ipt] + shift_nh;
      g_chisq_par1_nh_shift->SetPoint(ipt,valx,valy);
      cout << valy << " ";
    }
    cout << endl;
    for(int ipt=0;ipt<nbins1_ih;ipt++){      
      double valx = g_chisq_par1_ih->GetX()[ipt];
      double valy = g_chisq_par1_ih->GetY()[ipt] + shift_ih;
      g_chisq_par1_ih_shift->SetPoint(ipt,valx,valy);
      cout << valy << " ";
    }
    cout << endl;
    cout << "shift_nh = " << shift_nh << " - shift_ih = " << shift_ih << endl;
  }
  
  // Draw plots

  // Normal Hierarchy

  // 2-dimensional plots
  if(nDim==2){
    
    NuOscDmsqScheme_t orig_scheme = Opt::Instance()->gOscDmsqScheme; // save the scheme so we can put it back later
    
    if (!NuBarMode) {
      h_chisq_par2_par1_Ave_nh->GetXaxis()->SetTitle(util::ParamAsStringPaperRooTex(conv::kPhypMixing23).c_str());
      h_chisq_par2_par1_Ave_ih->GetXaxis()->SetTitle(util::ParamAsStringPaperRooTex(conv::kPhypMixing23).c_str());
      
      Opt::Instance()->gOscDmsqScheme = dmsqscheme_nh;
      h_chisq_par2_par1_Ave_nh->GetYaxis()->SetTitle(util::ParamAsStringPaperRooTex(conv::kPhypDm23sq).c_str());
      
      Opt::Instance()->gOscDmsqScheme = dmsqscheme_ih;
      h_chisq_par2_par1_Ave_ih->GetYaxis()->SetTitle(util::ParamAsStringPaperRooTex(conv::kPhypDm23sq).c_str());
    }
    else {
      h_chisq_par2_par1_Ave_nh->GetXaxis()->SetTitle(util::ParamAsStringPaperRooTex(conv::kPhypMixing23Bar).c_str());
      h_chisq_par2_par1_Ave_ih->GetXaxis()->SetTitle(util::ParamAsStringPaperRooTex(conv::kPhypMixing23Bar).c_str());
      
      Opt::Instance()->gOscDmsqScheme = dmsqscheme_nh;
      h_chisq_par2_par1_Ave_nh->GetYaxis()->SetTitle(util::ParamAsStringPaperRooTex(conv::kPhypDm23sqBar).c_str());
      
      Opt::Instance()->gOscDmsqScheme = dmsqscheme_ih;
      h_chisq_par2_par1_Ave_ih->GetYaxis()->SetTitle(util::ParamAsStringPaperRooTex(conv::kPhypDm23sqBar).c_str());
    }
    
    Opt::Instance()->gOscDmsqScheme = orig_scheme;
    
    TCanvas *c_chisq_par2_par1_nh = new TCanvas("c_chisq_par2_par1_nh","c_chisq_par2_par1_nh");
    h_chisq_par2_par1_Ave_nh->SetContour(100);
    h_chisq_par2_par1_Ave_nh->Draw("colz");
    gPad->Update();
    TPaletteAxis *palette_nh = (TPaletteAxis*)h_chisq_par2_par1_Ave_nh->GetListOfFunctions()->FindObject("palette");
    palette_nh->SetX1NDC(0.9);
    palette_nh->SetX2NDC(0.95);
    //palette_nh->SetY1NDC(0.);
    //palette_nh->SetY2NDC(.8);
    gPad->Modified();
    gPad->Update();
    //h_chisq_par2_par1_Ave_nh->DrawClone("colz");

    
    TCanvas *c_chisq_par2_par1_ih = new TCanvas("c_chisq_par2_par1_ih","c_chisq_par2_par1_ih");
    h_chisq_par2_par1_Ave_ih->SetContour(100);
    h_chisq_par2_par1_Ave_ih->Draw("colz");
    gPad->Update();
    TPaletteAxis *palette_ih = (TPaletteAxis*)h_chisq_par2_par1_Ave_ih->GetListOfFunctions()->FindObject("palette");
    palette_ih->SetX1NDC(0.9);
    palette_ih->SetX2NDC(0.95);
    // //palette_ih->SetY1NDC(0.);
    // //palette_ih->SetY2NDC(.8);
    gPad->Modified();
    gPad->Update();
    //h_chisq_par2_par1_Ave_ih->DrawClone("colz");

    
    TCanvas *c_cont_par2_par1 = new TCanvas("c_cont_par2_par1","c_cont_par2_par1");
    c_cont_par2_par1->cd();
    c_cont_par2_par1->SetGridx();
    c_cont_par2_par1->SetGridy();

    h_chisq_par2_par1_Ave_nh_Cont = (TH2F*)h_chisq_par2_par1_Ave_nh->Clone("");
    h_chisq_par2_par1_Ave_nh_Cont->SetContour(1);
    h_chisq_par2_par1_Ave_nh_Cont->SetContourLevel(0,CLval);
    h_chisq_par2_par1_Ave_nh_Cont->SetLineColor(kBlue);
    h_chisq_par2_par1_Ave_nh_Cont->Draw("cont3");
 
    h_chisq_par2_par1_Ave_ih_Cont = (TH2F*)h_chisq_par2_par1_Ave_ih->Clone("");
    h_chisq_par2_par1_Ave_ih_Cont->SetContour(1);
    h_chisq_par2_par1_Ave_ih_Cont->SetLineColor(kBlack);
    h_chisq_par2_par1_Ave_ih_Cont->SetContourLevel(0,CLval);
    h_chisq_par2_par1_Ave_ih_Cont->Draw("cont3 same");

    if(DrawLeg){
      TLegend *l_cont_par2_par1 = new TLegend(0.7,0.7,0.85,0.85);
      l_cont_par2_par1->SetFillColor(0);
      l_cont_par2_par1->SetLineColor(0);
      l_cont_par2_par1->AddEntry(h_chisq_par2_par1_Ave_nh,LegendName_nh,"L");
      l_cont_par2_par1->AddEntry(h_chisq_par2_par1_Ave_ih,LegendName_ih,"L");
      l_cont_par2_par1->Draw();
    }
  }
  // 1-dimensional plots
  else if(nDim==1){
    
    TCanvas *c_chisq_par1_nh = new TCanvas("c_chisq_par1_nh","c_chisq_par1_nh");
    double min = TMath::Min(0.,(double)g_chisq_par1_nh_Err->GetMinimum());
    min = TMath::Max(0.,min);
    double max = g_chisq_par1_nh_Err->GetMaximum() * 1.25;
    g_chisq_par1_nh_Err->SetMinimum(max);
    g_chisq_par1_nh_Err->SetMinimum(0);
    g_chisq_par1_nh_Err->SetFillColor(kRed);
    g_chisq_par1_nh_Err->SetFillStyle(3001);
    g_chisq_par1_nh_Err->DrawClone("a2");
    g_chisq_par1_nh->SetFillColor(kBlue);
    g_chisq_par1_nh->DrawClone("c same");

    TCanvas *c_chisq_par1_ih = new TCanvas("c_chisq_par1_ih","c_chisq_par1_ih");
    double min = TMath::Min(0.,(double)g_chisq_par1_ih_Err->GetMinimum());
    min = TMath::Max(0.,min);
    double max = g_chisq_par1_ih_Err->GetMaximum() * 1.25;
    g_chisq_par1_ih_Err->SetMinimum(max);
    g_chisq_par1_ih_Err->SetMinimum(0);
    g_chisq_par1_ih_Err->SetFillColor(kRed);
    g_chisq_par1_ih_Err->SetFillStyle(3001);
    g_chisq_par1_ih_Err->DrawClone("a2");
    g_chisq_par1_ih->SetFillColor(kBlue);
    g_chisq_par1_ih->DrawClone("c same");

    if(doShiftMin){
      TCanvas *c_chisq_par1_mhall = new TCanvas("c_chisq_par1_mhall","c_chisq_par1_mhall");
      g_chisq_par1_nh_shift->SetMaximum(10);
      g_chisq_par1_nh_shift->SetMinimum(0);
      g_chisq_par1_nh_shift->SetLineColor(kBlack);
      g_chisq_par1_nh_shift->SetMarkerColor(kBlack);
      g_chisq_par1_nh_shift->SetLineWidth(2);
      if(*parname1_nh=="deltacp") g_chisq_par1_nh_shift->GetXaxis()->SetTitle("#delta_{CP}");
      g_chisq_par1_nh_shift->GetYaxis()->SetTitle("#Delta #chi^{2}");

      g_chisq_par1_nh_shift->DrawClone("ac");
      //g_chisq_par1_nh_shift->DrawClone("ap");
      g_chisq_par1_ih_shift->SetLineColor(kRed);
      g_chisq_par1_ih_shift->SetMarkerColor(kRed);
      g_chisq_par1_ih_shift->SetLineWidth(2);
      g_chisq_par1_ih_shift->GetXaxis()->SetTitle(parname2_nh[0].c_str());
      g_chisq_par1_ih_shift->DrawClone("c same");
      //g_chisq_par1_ih_shift->DrawClone("p same");
      
      if(DrawLeg){
        TLegend *l_cont_par1 = new TLegend(0.7,0.7,0.85,0.85);
        l_cont_par1->SetFillColor(0);
        l_cont_par1->SetLineColor(0);
        l_cont_par1->AddEntry(g_chisq_par1_nh_shift,LegendName_nh,"L");
        l_cont_par1->AddEntry(g_chisq_par1_ih_shift,LegendName_ih,"L");
        l_cont_par1->Draw();
      }
    }
  
  } // NDim == 1


  // Print plots
  
  if(PrintPlots){
    
    char name[200];
    
    if(nDim==1){
    }
    else if(nDim==2){
      sprintf(name,"%s%s%s","Like2D_NH_" ,tag,img_format.c_str());
      c_chisq_par2_par1_nh->Print(name);
      sprintf(name,"%s%s%s","Like2D_IH_",tag,img_format.c_str());
      c_chisq_par2_par1_ih->Print(name);
      sprintf(name,"%s%s%s","Cont2D_"   ,tag,img_format.c_str());
      c_cont_par2_par1->Print(name);
    }
    
  }

  // Write histograms to root file
  if (output_file) {
    TFile hist_file(output_file,"UPDATE");
    hist_file.cd();
    if(nDim==2){
      h_chisq_par2_par1_Ave_nh->Write();
      h_chisq_par2_par1_Ave_ih->Write();
    }
    if(nDim==1) {
      h_chisq_par1_Ave_nh->Write();
      h_chisq_par1_Ave_ih->Write();
    }
    
    cout<<"Histograms written to "<<output_file<<endl;
  }
      

  return;
}

  
void set_contour_style(TH2D *contour, Int_t color, Int_t width, Int_t style){
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
  contour->SetContour(1);//number of contour levels                                                                                                                                    
  contour->SetContourLevel(0,0.0);//SetContourLevel(int level, double value)                                                                                                           
  contour->SetStats(0);
} 

TGraph * get_contour(TH2D * h) {
  TCanvas * c = new TCanvas();
  h->Draw("cont list");
  c->Update();
  TObjArray * contours  = gROOT->GetListOfSpecials()->FindObject("contours");
  if(contours->GetSize() != 1)
    exit(1);
  TList     * list      = (TList*)contours->At(0);
  if(list->GetSize() != 1)
    exit(1);
  TGraph * graph = (TGraph*)list->First();
  delete c;
  return graph;
}
