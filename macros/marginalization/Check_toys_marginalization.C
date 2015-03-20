
//
// ------------------------------------------------------------------------------------------------
//
// Check the extraction of the systematic and oscillation parameters done for marginalization
// as well as the Ereco distributions.
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------



using namespace std;
using namespace valor;
using namespace valor::conv;

void Check_toys_marginalization
(
 const char * outfilename          = "ChisqMarg",
 
 const int    IDoutfile            = 0,
 const int    first_bin_freepar0   = 0,   // first freepar of the sub-range
 const int    last_bin_freepar0    = 4, // last freepar of the sub-range
 const int    first_bin_freepar1   = 0,   // first freepar of the sub-range
 const int    last_bin_freepar1    = 4, // last freepar of the sub-range
 
 //const char * infilename           = "../../inputs/marginalization/toys/mixing23_1p7E20_postfit_NoRethrow_Merged",
 //const char * infilename           = "MH_FitMarg_mixing23_1p7E20_prefit",
 //const char * infilename           = "MH_FitMarg_mixing23_1p7E20_prefit_Merged",
 //const char * infilename           = "submit_mixing23_1p7E20_postfit_NoRethrow",
 const char * infilename           = "submit_mixing23_1p7E20_postfit_NoRethrow_Merged",
 //const char * infilename           = "../../inputs/marginalization/toys/mixing23_dmsq32_1p7E20_postfit_NoRethrow_Merged",
 
 const char * config_file          = "analysis_configuration_1.7E20.xml",
 //const char * config_file          = "analysis_configuration_1.7E20_prefit.xml",   /* analysis configuration file */
 //const char * config_file          = "analysis_configuration_1.7E20.xml",   /* analysis configuration file */
 
 //const char * fit_param_list       = "nuosc3fbar_fixnonbar_fix5+nuisance",
 const char * fit_param_list       = "nuosc3fbar_fixnonbar_12fix+nuisance",
 //const char * fit_param_list      = "nuosc3fbar_numudisapp_constrainall+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",
 
 const int    NtoysMarg            = 1000,   // # toys marginalization
 
 const char * freeparname0        = "mixing23Bar",  // free parameter
 const double freeparmin0         = 0,              // minimum value of free parameter
 const double freeparmax0         = 1,            // maxinum value of free parameter
 const int    NptScan0            = 21,            // # pts of free parameter to scan
 
 //const char * freeparname1        = "",  // free parameter
 const char * freeparname1        = "dmsq23Bar",  // free parameter
 const double freeparmin1         = 0.0002, // minimum value of free parameter
 const double freeparmax1         = 0.02, // maxinum value of free parameter
 const int    NptScan1            = 21,          // # pts of free parameter to scan
 
 int         tweak_systematics    = 1,       /* generate toy MC with randomized systematics? */
 
 const bool  doMargMH             = true,    // marginalize Mass Hierarchy
 const bool  doRethrowUnphysical  = false,    // rethrow unphysical systematic toys
 
 const char * dmsqscheme          = "abs23",      /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy           = "inverted",     /* normal or inverted */
 const char * mixparamchoice      = "sinsqtheta", /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant      = "one",       /* one==kNuOctantOne two==kNuOctantTwo */
 
 const char * xmlfile = "fakedata/201411.fakedata/fakedata_201411_0p170e21_A.xml",  // input XML file
 bool using_xml_file  = true,
 bool draw_toys       = true
 )
{
    
    FileStat_t dummy;
    
    TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);
    
    // Check for config file, and XML data file.
    // If not found in inputs or config, see if it works as a normal path
    if (gSystem->GetPathInfo(config.Data(), dummy))
    {
        config = config_file;
        if (gSystem->GetPathInfo(config.Data(), dummy)) {
            cerr<<"Config file not found: "<<config_file<<endl;
            exit(-1);
        }
    }
    
    if (gSystem->GetPathInfo(xmlfile2.Data(), dummy))
    {
        xmlfile2 = xmlfile;
        if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
            cerr<<"Data input file not found: "<<xmlfile<<endl;
            exit(-1);
        }
    }
    
    int ifakedata = 0;
    SuperkObservation * skobs_fakedata = new SuperkObservation;
    skobs_fakedata->Read(xmlfile2.Data(),ifakedata,using_xml_file);
    
    // Global options
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
    Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
    Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
    Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
    if(Opt::Instance()->Undefined()) {
        cerr << "Global options not defined correctly. Exiting..." << endl;
        cerr << *Opt::Instance() << endl;
        exit(1);
    }
    Opt::Instance()->SetDefaultParameterValues();
    
    // Get superk spectrum predictor & configure
    SuperkExpectation * skexp = new SuperkExpectation;
    skexp->Read(config.Data());
    
    
    PhysParamId_t freepar0;
    string sfreeparname0(freeparname0);
    // NuBar mode
    if(sfreeparname0=="mixing23Bar")      freepar0 = conv::kPhypMixing23Bar;
    else if(sfreeparname0=="mixing13Bar") freepar0 = conv::kPhypMixing13Bar;
    else if(sfreeparname0=="dmsq23Bar")   freepar0 = conv::kPhypDm23sqBar;
    else if(sfreeparname0=="dmsq12Bar")   freepar0 = conv::kPhypDm12sqBar;
    else if(sfreeparname0=="mixing12Bar") freepar0 = conv::kPhypMixing12Bar;
    else if(sfreeparname0=="deltacpBar")  freepar0 = conv::kPhypCPVBar;
    // Nu mode
    else if(sfreeparname0=="mixing23")    freepar0 = conv::kPhypMixing23;
    else if(sfreeparname0=="mixing13")    freepar0 = conv::kPhypMixing13;
    else if(sfreeparname0=="dmsq23")      freepar0 = conv::kPhypDm23sq;
    else if(sfreeparname0=="dmsq12")      freepar0 = conv::kPhypDm12sq;
    else if(sfreeparname0=="mixing12")    freepar0 = conv::kPhypMixing12;
    else if(sfreeparname0=="deltacp")     freepar0 = conv::kPhypCPV;
    //
    else{
        freepar0 = conv::kPhypUndefined;
        cout << "Param 0 is not defined" << endl;
        exit(1);
    }
    
    PhysParamId_t freepar1 = conv::kPhypUndefined;
    string sfreeparname1(freeparname1);
    // NuBar mode
    if(sfreeparname1=="mixing23Bar")      freepar1 = conv::kPhypMixing23Bar;
    else if(sfreeparname1=="mixing13Bar") freepar1 = conv::kPhypMixing13Bar;
    else if(sfreeparname1=="dmsq23Bar")   freepar1 = conv::kPhypDm23sqBar;
    else if(sfreeparname1=="dmsq12Bar")   freepar1 = conv::kPhypDm12sqBar;
    else if(sfreeparname1=="mixing12Bar") freepar1 = conv::kPhypMixing12Bar;
    else if(sfreeparname1=="deltacpBar")  freepar1 = conv::kPhypCPVBar;
    // Nu mode
    else if(sfreeparname1=="mixing23")    freepar1 = conv::kPhypMixing23;
    else if(sfreeparname1=="mixing13")    freepar1 = conv::kPhypMixing13;
    else if(sfreeparname1=="dmsq23")      freepar1 = conv::kPhypDm23sq;
    else if(sfreeparname1=="dmsq12")      freepar1 = conv::kPhypDm12sq;
    else if(sfreeparname1=="mixing12")    freepar1 = conv::kPhypMixing12;
    else if(sfreeparname1=="deltacp")     freepar1 = conv::kPhypCPV;
    //
    else{
        freepar1 = conv::kPhypUndefined;
        cout << "1D fit" << endl;
    }
    
    
    // Set free param, ntoys, # pt to scan, min par, max par (if min and max not passed, use default ones for that par)
    
    // free parameter
    PhysParamId_t parfree[2];
    parfree[0] = freepar0;
    parfree[1] = freepar1;
    
    // # of toys for marginalization
    //int ntoys = NtoysMarg;
    int ntoys = 1;
    
    // # of free parameter pts to scan
    int nbinspar[2];
    nbinspar[0] = NptScan0;
    nbinspar[1] = NptScan1;
    
    // free parameter range
    double minpar[2];
    double maxpar[2];
    minpar[0] = freeparmin0;
    maxpar[0] = freeparmax0;
    minpar[1] = freeparmin1;
    maxpar[1] = freeparmax1;
    
    // free parameter sub-range
    int subminbin[2];
    int submaxbin[2];
    subminbin[0] = first_bin_freepar0; submaxbin[0] = last_bin_freepar0;
    subminbin[1] = first_bin_freepar1; submaxbin[1] = last_bin_freepar1;
    
    
    // Implement the constructor giving the following input parameters:
    FitMarginalization *fitmarg = new FitMarginalization();
    fitmarg->SetPhysParFree(parfree,NtoysMarg,nbinspar,minpar,maxpar);
    fitmarg->SetSubRange(subminbin,submaxbin); // By default it scans all the range
    fitmarg->Init(); // Before need SetPhysParFree() to initialize the arrays
    
    // Set fit parameter list
    fitmarg->SetFitParamList(fit_param_list);
    fitmarg->PrintPhysHypothesis();
    
    // Set Mass Hierarchy marginalization
    fitmarg->SetMargMassHierarchy(doMargMH); // false by default
    
    // Decide whether to rethrow or not the unphysical systematic toys
    fitmarg->SetRethrowSystUnphysical(doRethrowUnphysical);
    
    // Save parameter throws
    fitmarg->SetSaveParamLists();
    
    // Produce toys and fill marginalized chisq
    
    // Set priors for not fixed oscillation parameters
    fitmarg->SetPrior(kPhypMixing23Bar,"flat",0,1);
    fitmarg->SetPrior(kPhypMixing13Bar,"flat",0,1);
    fitmarg->SetPrior(kPhypMixing12Bar,"flat",0,1);
    fitmarg->SetPrior(kPhypDm23sqBar,"flat",conv::kMinFitDm23sq,conv::kMaxFitDm23sq);
    fitmarg->SetPrior(kPhypDm12sqBar,"flat",conv::kMinFitDm12sq,conv::kMaxFitDm12sq);
    fitmarg->SetPrior(kPhypCPVBar,"flat",-TMath::Pi(),+TMath::Pi());
    
    
    // In order to get toy extraction you need to
    //for(int itoy=0;itoy<NtoysMarg;itoy++){
    
    fitmarg->MakeToysMarginalization(tweak_systematics,skexp);
    
    
    //}
    
    //  fitmarg->GetToyMargHisto(conv::kSKsRHCFCFV1RingMuLike,0,0,0,0)->Draw();
    //fitmarg->GetToyMargHisto(conv::kSKsRHCFCFV1RingELike,0,0,0,0)->Draw();
    
    ///
    if(draw_toys)
    {
        valor::conv::SuperkSample_t s = conv::kSKsRHCFCFV1RingMuLike; ///Set the sample
        
        int ipar0 = first_bin_freepar0; //set which bin to draw toys for
        int ipar1 = first_bin_freepar1; // set which bin to draw toys for
        
        double step0 = (freeparmax0 - freeparmin0)/NptScan0;
        double step1 = (freeparmax1 - freeparmin1)/NptScan1;
        
        double val0 = freeparmin0 + ipar0 * step0; // freeparmin0 + number0 * ipar0;
        double val1 = freeparmin1 + ipar1 * step1; // freeparmin1 + number1 * ipar1;
        
        TH1F* h_marg_toy;
        TH1F* h_up = new TH1F("h_up","h_up",util::NumErecoBins(s),util::ErecoBinEdges(s));
        TH1F* h_down = new TH1F("h_down","h_down",util::NumErecoBins(s),util::ErecoBinEdges(s));
        TH1D* h_nominal;
        TH2F* h_2d = new TH2F("h_2d","Envelope",util::NumErecoBins(s),util::ErecoBinEdges(s),1000,0,6);
        
        //Normal Hierarchy
        if(doMargMH || Opt::Instance()->gMassHierarchyAssumed == kNmhNormal){
        for(int toy=0;toy<NtoysMarg;toy++){
            h_marg_toy = fitmarg->GetToyMargHisto(s,toy,ipar0-subminbin[0],ipar1-subminbin[1],0);
            for(int i=1;i<=util::NumErecoBins(s);i++)
            {
                h_2d->Fill(h_marg_toy->GetBinCenter(i),h_marg_toy->GetBinContent(i));
            }
        }
        }
        
        //Inverted Hierarchy
        if(doMargMH || Opt::Instance()->gMassHierarchyAssumed == kNmhInverted){
        for(int toy=0;toy<NtoysMarg;toy++){
            h_marg_toy = fitmarg->GetToyMargHisto(s,toy,0,0,1);
            for(int i=1;i<=util::NumErecoBins(s);i++)
            {
                h_2d->Fill(h_marg_toy->GetBinCenter(i),h_marg_toy->GetBinContent(i));
            }
        }
        }
        
        for(int i=1;i<=util::NumErecoBins(s);i++)
        {
            TH1D* h_proj = h_2d->ProjectionY("",i,i);
            double average = average(h_proj);
            double rms = rms(h_proj,average);
            
            h_up->SetBinContent(i,average+ rms);
            h_down->SetBinContent(i,average-rms);
            
        }
        
        skexp->ResetSystematicVariations();
        FitParamList * fparm  = util::NamedParameterList(fit_param_list); //only used to get nominal
        
        fparm->Update(conv::kPhypMixing12,        conv::kGlobFitSinsqTheta12);
        fparm->Update(conv::kPhypMixing13,        conv::kGlobFitSinsqTheta13);
        fparm->Update(conv::kPhypMixing23,        conv::kGlobFitSinsqTheta23);
        fparm->Update(conv::kPhypCPV    ,         conv::kGlobFitDeltaCP);
        fparm->Update(conv::kPhypDm12sq,          conv::kGlobFitDm12sq);
        fparm->Update(conv::kPhypDm23sq,          conv::kGlobFitDm23sq);
        fparm->Update(conv::kPhypMixing12Bar,     conv::kGlobFitSinsqTheta12);
        fparm->Update(conv::kPhypMixing13Bar,     conv::kGlobFitSinsqTheta13);
        fparm->Update(conv::kPhypMixing23Bar,     conv::kGlobFitSinsqTheta23);
        fparm->Update(conv::kPhypCPVBar    ,      conv::kGlobFitDeltaCP);
        fparm->Update(conv::kPhypDm12sqBar,       conv::kGlobFitDm12sq);
        fparm->Update(conv::kPhypDm23sqBar,       conv::kGlobFitDm23sq);
        fparm->Update(freepar0,val0);
        fparm->Update(freepar1,val1);
        
        cout << util::AsString(freepar0) << ":\t" << val0 << util::AsString(freepar1) << ":\t" << val1 << endl;
        
        skexp->Generate(PhysHypothesis(*fparm));
        h_nominal = skexp->CurrPredictedErecoSpectrum(s);
        
        TCanvas * c = new TCanvas("c","",20,20,1000,1000);
        h_2d->Draw("colz");
        h_up->Draw("same");
        h_down->Draw("same");
        h_nominal->Draw("same");
        
        c->SaveAs("draw_marg_toys.root");
        
    }
    
    // Store toy distributions in output file  
    char outname[200];
    sprintf(outname,"check_toy_throws_%d.root",outfilename,IDoutfile);
    //~ cout << outname << endl;
    TFile *fout = new TFile(outname,"RECREATE");
    fout->cd();
    TObjArray * arr = fitmarg->GetToyParamThrowHistograms();
    arr->Write();
    fout->Close();
  
    
    // // Store ChisqMarg (full or partial) in output file
    // char outname[200];
    // sprintf(outname,"%s_%d.root",outfilename,IDoutfile);
    // cout << outname << endl;
    // TFile *fout = new TFile(outname,"RECREATE");
    // fitmarg->GetChisqMarg()->Write();
    // fitmarg->GetChisqMarg2D()->Write();
    // fout->Close();
    
    return;
}

double average(TH1D* h1)
{
    double avg(0);
    
    for(int i=1; i<=h1->GetNbinsX();i++)
    {
        avg += h1->GetBinCenter(i)*h1->GetBinContent(i);
    }
    
    avg = avg/h1->GetEntries();
    
    return avg;
}

double rms(TH1D* h1,double avg)
{
    double rms(0);
    
    for(int i=1; i<=h1->GetNbinsX();i++)
    {
        rms += ((h1->GetBinCenter(i) - avg) **2) * (h1->GetBinContent(i));
    }
    
    rms = rms/h1->GetEntries();
    rms = rms**0.5;
    return rms;
}
