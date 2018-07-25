#include <vector>
#include <iostream>
#include <algorithm>
#include <utility>

using namespace std;
using namespace valor;
using namespace valor::conv;

const double CLfrac_1sigma = 0.6827;  // 1sigma% CL
const double CLfrac_90CL   = 0.9;  // 90% CL
const double CLfrac_95CL   = 0.95;  // 95% CL
const double CLfrac_2sigma = 0.9545;  // 2sigma% CL

double chi2(double* xx, double* p)
{
   double x = xx[0];
   double dof2 = 0.5 * p[0];
   double numerator = TMath::Exp(-0.5 * x) * TMath::Power(x, dof2 - 1.0);
   double denominator = TMath::Power(2.0, dof2) * TMath::Gamma(dof2);
   return numerator / denominator;
}

void plot_brazilian_mixing23
(
   const string tag = "dcp_run18_react_Data",
   TString titleX = "#delta_{CP} (rad)",
   TString titleY = "-2#Deltaln(L)",
   const bool doGlobalMin  = true,
   const bool draw_preliminary = true,
   TString str_preliminary = "T2K Run1-9c Preliminary",
   const double truepar0 = -1.833,
   const string truemh0 = "normal",
   const double exclusion_value = 0.0,
   const char *infilename_toys = "/home/phrdqd/VaLOR/analysis/tn/final/fitmarg_fc_dcp_normal_final_10k_-1.833_FullChisq.root",
   const char *fit_file_nh = "/home/phrdqd/VaLOR/analysis/tn/final/fit_oa2017_data_normal_reactor_deltacp_nexp_1_-1.601_-999.000-job0-_0_0.root",
   const char *fit_file_ih = "/home/phrdqd/VaLOR/analysis/tn/final/fit_oa2017_data_inverted_reactor_deltacp_nexp_1_-1.601_-999.000-job0-_0_0.root",
   const char *fcname = "/home/phrdqd/VaLOR/analysis/tn/final/FCgraphs.root"
)
{
   str_preliminary = "T2K Run1-9c Preliminary";
   gStyle->SetOptTitle(0);
   const int Nfiles = 1;
   const int Nmh = 2;
   char *infilename[Nfiles][Nmh] =
   {
      {
         fit_file_nh,
         fit_file_ih
      }
   };
    
   char *latexname[Nfiles][Nmh] =
   {
      {"Normal","Inverted"}
   };

  char name[200];
  
  TFile *filein[Nfiles][Nmh];
  TTree *treein[Nfiles][Nmh];
  
  TH1D *hchisq1D[Nfiles][Nmh];
  TH1D *hdeltachisq1D[Nfiles][Nmh];
  TH1D *hdeltachisq1D_GlobMin[Nfiles][Nmh];

  double bestfit_x = -999.;
  double bestfit_mh = -999.;
  
  TCanvas *cChisq; //[Nfiles][Nmh];
  
  for(int ifile=0;ifile<Nfiles;ifile++)
  {
     for(int imh=0;imh<Nmh;imh++)
     {
        cout << "Get file " << infilename[ifile][imh] << endl;
        filein[ifile][imh] = new TFile(infilename[ifile][imh],"READ");
        if(!filein[ifile][imh]->IsOpen())
        {
           cerr << "File " << infilename[ifile][imh] << " is not open!" << endl;
           exit(1);
        }
        filein[ifile][imh]->cd();
        treein[ifile][imh] = (TTree*)filein[ifile][imh]->Get("fitnt");
      
        sprintf(name,"hchisq1D_file%d_mh%d",ifile,imh);
        hchisq1D[ifile][imh] = new TH1D(name,name,1,0,1);
        treein[ifile][imh]->SetBranchAddress("hchisq1D", &hchisq1D[ifile][imh] );
      
        sprintf(name,"hdeltachisq1D_file%d_mh%d",ifile,imh);
        hdeltachisq1D[ifile][imh] = new TH1D(name,name,1,0,1);
        treein[ifile][imh]->SetBranchAddress("hdeltachisq1D", &hdeltachisq1D[ifile][imh] );
      
        // Get event
        treein[ifile][imh]->GetEntry(0);
      
        // Print chisq distributions
        cout << std::setprecision(9);
        cout << "Chisq file " << ifile << " - mh = " << imh << " :" << endl;
        cout << "# of entries = " << hchisq1D[ifile][imh]->GetEntries() << endl;
        for(int i=0;i<hchisq1D[ifile][imh]->GetNbinsX();i++)
        {
           cout << hchisq1D[ifile][imh]->GetBinContent(i+1) << " ";
        }
        cout << endl;
        cout << "DeltaChisq: # of entries = " << hdeltachisq1D[ifile][imh]->GetEntries() << endl;
        for(int i=0;i<hdeltachisq1D[ifile][imh]->GetNbinsX();i++)
        {
           cout << hdeltachisq1D[ifile][imh]->GetBinContent(i+1) << " ";
        }
        cout << endl;
      
        sprintf(name,"hdeltachisq1D_GlobMin_file%d_mh%d",ifile,imh);
        hdeltachisq1D_GlobMin[ifile][imh] = (TH1D*)hchisq1D[ifile][imh]->Clone(name);
        hdeltachisq1D_GlobMin[ifile][imh]->Reset();
      } // end loop Nmh
    
      // Shift the chisq wrt to the global minimum
      double ChisqMin_0 = hchisq1D[ifile][0]->GetBinContent(hchisq1D[ifile][0]->GetMinimumBin());
      double ChisqMin_1 = hchisq1D[ifile][1]->GetBinContent(hchisq1D[ifile][1]->GetMinimumBin());
      double ChisqMin_Glob = ChisqMin_0;
      if( ChisqMin_1 < ChisqMin_Glob )
      {
         ChisqMin_Glob = ChisqMin_1;
      }
    
      int XminBin_0 = hchisq1D[ifile][0]->GetMinimumBin();
      double Xmin_0 = hchisq1D[ifile][0]->GetBinCenter(XminBin_0); 
      int XminBin_1 = -999;
      double Xmin_1 = -999;
      if(hchisq1D[ifile][1])
      {
         XminBin_1 = hchisq1D[ifile][1]->GetMinimumBin();
         Xmin_1 = hchisq1D[ifile][1]->GetBinCenter(XminBin_1);
      }
    
      if(ChisqMin_0<ChisqMin_1)
      {
         cout << "Best-fit = " << Xmin_0 << endl;
      }
      else
      {
         cout << "Best-fit = " << Xmin_1 << endl;
      }
      cout << "ChisqMin Glob = " << ChisqMin_Glob << endl;
    
      for(int imh=0;imh<Nmh;imh++)
      {
         double ChisqMin_Curr = hchisq1D[ifile][imh]->GetBinContent(hchisq1D[ifile][imh]->GetMinimumBin());
         for(int i=0;i<hdeltachisq1D_GlobMin[ifile][imh]->GetNbinsX();i++)
         {
            double cont = hchisq1D[ifile][imh]->GetBinContent(i+1);
            if(doGlobalMin)
            {
               hdeltachisq1D_GlobMin[ifile][imh]->SetBinContent(i+1,cont-ChisqMin_Glob);
            }
            else
            {
               hdeltachisq1D_GlobMin[ifile][imh]->SetBinContent(i+1,cont-ChisqMin_Curr);
            }
         }
      }
   } // end loop Nfiles
  
   // Get critical values
   TGraph *gFC1sigma_nh;
   TGraph *gFC1sigma_ih;
   TGraph *gFC90CL_nh;
   TGraph *gFC90CL_ih;
   TGraph *gFC2sigma_nh;
   TGraph *gFC2sigma_ih;
   TGraph *gFC3sigma_nh;
   TGraph *gFC3sigma_ih;
   TFile *fin_fc_nh;
   TFile *fin_fc_ih;
  
   fin_fc = new TFile(fcname,"READ");
    
   gFC90CL_nh      = (TGraph*) fin_fc->Get("gFC90CL_nh");
   gFC2sigma_nh    = (TGraph*) fin_fc->Get("gFC2sigma_nh");
   gFC3sigma_nh    = (TGraph*) fin_fc->Get("gFC3sigma_nh");
   gFC1sigma_nh    = (TGraph*) fin_fc->Get("gFC1sigma_nh");

   gFC90CL_ih      = (TGraph*) fin_fc->Get("gFC90CL_ih");
   gFC2sigma_ih    = (TGraph*) fin_fc->Get("gFC2sigma_ih");
   gFC3sigma_ih    = (TGraph*) fin_fc->Get("gFC3sigma_ih");
   gFC1sigma_ih    = (TGraph*) fin_fc->Get("gFC1sigma_ih");

   TGaxis::SetMaxDigits(2);
    
   TFile *filetoys = new TFile(infilename_toys,"READ");
   if(!filetoys->IsOpen())
   {
      cerr << "File " << infilename << " is not open" << endl;
   }
   cout << "Reading file " << infilename_toys << endl;
   TTree *treetoys = (TTree*)filetoys->Get("fitnt");
   if(!treetoys)
   {
      cerr << "Tree object is empty" << endl;
   }
    
   TH1D *hdeltachisq1D_toy_nh; // take it from hchisq1D_toy_nh at each event
   TH1D *hdeltachisq1D_toy_ih; // take it from hchisq1D_toy_nh at each event
    
   TH1D *hNevts_1Rmu    = new TH1D("hNevts_1Rmu","hNevts_1Rmu",500,50,200);
   TH1D *hNevts_RHC1Rmu = new TH1D("hNevts_RHC1Rmu","hNevts_RHC1Rmu",240,20,100);
   TH1D *hNevts_1Re    = new TH1D("hNevts_1Re","hNevts_1Re",120,0,60);
   TH1D *hNevts_RHC1Re = new TH1D("hNevts_RHC1Re","hNevts_RHC1Re",40,0,20);
    
   TH1D *hchisq1D_toy_nh = new TH1D("hchisq1D_toy_nh","hchisq1D_toy_nh",1,0,1);
   TH1D *hchisq1D_toy_ih = new TH1D("hchisq1D_toy_ih","hchisq1D_toy_ih",1,0,1);
   TH2D *hErecoToy_1Rmu = new TH2D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
   TH2D *hErecoToy_1Re = new TH2D("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
   TH2D *hErecoToy_RHC1Rmu = new TH2D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
   TH2D *hErecoToy_RHC1Re = new TH2D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);
    
   treetoys->SetBranchAddress("hchisq1D_nh",        &hchisq1D_toy_nh      );
   treetoys->SetBranchAddress("hchisq1D_ih",        &hchisq1D_toy_ih      );
   treetoys->SetBranchAddress("hErecoToy_1Rmu",     &hErecoToy_1Rmu     );
   treetoys->SetBranchAddress("hErecoToy_1Re",      &hErecoToy_1Re      );
   treetoys->SetBranchAddress("hErecoToy_RHC1Rmu",  &hErecoToy_RHC1Rmu  );
   treetoys->SetBranchAddress("hErecoToy_RHC1Re",   &hErecoToy_RHC1Re   );
    
   int Nentries = treetoys->GetEntries();
   treetoys->GetEntry(0);
   int NbinsChisq = hchisq1D_toy_nh->GetNbinsX();
    
   vector<double> *vec_deltachisq_bin_nh = new vector<double>[NbinsChisq];
   vector<double> *vec_deltachisq_bin_ih = new vector<double>[NbinsChisq];
   vector<double> *vec_like_penalty      = new vector<double>[NbinsChisq];
    
   TH1D **hchisq1D_toy_nh_bin      = new TH1D[NbinsChisq];
   TH1D **hchisq1D_toy_ih_bin      = new TH1D[NbinsChisq];
   TH1D **hdeltachisq1D_toy_nh_bin = new TH1D[NbinsChisq];
   TH1D **hdeltachisq1D_toy_ih_bin = new TH1D[NbinsChisq];
   for(int ibin=0;ibin<NbinsChisq;ibin++)
   {
      hchisq1D_toy_nh_bin[ibin] = new TH1D(Form("hchisq1D_toy_nh_bin%d",ibin),
            Form("hchisq1D_toy_nh_bin%d",ibin),500000,0,500);
      hchisq1D_toy_ih_bin[ibin] = new TH1D(Form("hchisq1D_toy_ih_bin%d",ibin),
            Form("hchisq1D_toy_ih_bin%d",ibin),500000,0,500);
      hdeltachisq1D_toy_nh_bin[ibin] = new TH1D(Form("hdeltachisq1D_toy_nh_bin%d",ibin),
            Form("hdeltachisq1D_toy_nh_bin%d",ibin),500000,0,35);
      hdeltachisq1D_toy_ih_bin[ibin] = new TH1D(Form("hdeltachisq1D_toy_ih_bin%d",ibin),
            Form("hdeltachisq1D_toy_ih_bin%d",ibin),500000,0,35);
   }
    
   // mean sensitivity
   TH1D *hchisq1D_nh_mean = hchisq1D_toy_nh->Clone("hchisq1D_nh_mean");
   hchisq1D_nh_mean->Reset();
   TH1D *hchisq1D_ih_mean = hchisq1D_toy_ih->Clone("hchisq1D_ih_mean");
   hchisq1D_ih_mean->Reset();
    
   int nbins = hchisq1D_ih_mean->GetXaxis()->GetNbins();
   double lowedge = hchisq1D_ih_mean->GetXaxis()->GetBinLowEdge(1);
   double upedge = hchisq1D_ih_mean->GetXaxis()->GetBinUpEdge(nbins);
   TH2D *hdcp_Vs_deltachisq_nh = new TH2D("hdcp_Vs_deltachisq_nh","hdcp_Vs_deltachisq_nh",
         nbins,lowedge,upedge,100,0,40);
   TH2D *hdcp_Vs_deltachisq_ih = new TH2D("hdcp_Vs_deltachisq_ih","hdcp_Vs_deltachisq_ih",
         nbins,lowedge,upedge,100,0,40);

   cout << "Enter the loop of " << Nentries << " entries" << endl;
    
   for(int ientry=0;ientry<Nentries;ientry++)
   {
      treetoys->GetEntry(ientry);
      if(ientry%1000==0)
      {
         cout << "Entry " << ientry << endl;
      }
      
      // Take global minimum
      double min_glob = 0.;
      int min_hierarchy = 0;
      double min_dcp = -999;
      double min_nh = hchisq1D_toy_nh->GetBinContent(hchisq1D_toy_nh->GetMinimumBin());
      double min_ih = hchisq1D_toy_ih->GetBinContent(hchisq1D_toy_ih->GetMinimumBin());
      if(min_nh < min_ih)
      {
         min_glob = min_nh;
         min_hierarchy = 1;
         min_dcp = hchisq1D_toy_nh->GetBinCenter(hchisq1D_toy_nh->GetMinimumBin());
      }
      else
      {
         min_glob = min_ih;
         min_hierarchy = 2;
         min_dcp = hchisq1D_toy_nh->GetBinCenter(hchisq1D_toy_ih->GetMinimumBin());
      }

      // Build mean sensitivity
      hchisq1D_nh_mean->Add(hchisq1D_toy_nh,1.);
      hchisq1D_ih_mean->Add(hchisq1D_toy_ih,1.);
      
      // Take global deltachisq and shift deltachisq of nh and ih
      hdeltachisq1D_toy_nh = (TH1D*) hchisq1D_toy_nh ->Clone("hchisq1D_toy_nh");
      hdeltachisq1D_toy_ih = (TH1D*) hchisq1D_toy_ih ->Clone("hchisq1D_toy_ih");
      
      for(int ibin=0;ibin<NbinsChisq;ibin++)
      {  
         double bincont = hdeltachisq1D_toy_nh->GetBinContent(ibin+1);
         hdeltachisq1D_toy_nh->SetBinContent(ibin+1,bincont-min_glob);
         bincont = hdeltachisq1D_toy_ih->GetBinContent(ibin+1);
         hdeltachisq1D_toy_ih->SetBinContent(ibin+1,bincont-min_glob);
      }
      
      // Fill deltaCP bin histograms
      
      for(int ibin=0;ibin<NbinsChisq;ibin++)
      {
         double chisqval = hchisq1D_toy_nh->GetBinContent(ibin+1);
         hchisq1D_toy_nh_bin[ibin]->Fill(chisqval);
         double deltachisqval = hdeltachisq1D_toy_nh->GetBinContent(ibin+1);
         hdeltachisq1D_toy_nh_bin[ibin]->Fill(deltachisqval);
         vec_deltachisq_bin_nh[ibin].push_back(deltachisqval);
        
         hdcp_Vs_deltachisq_nh->Fill(hdeltachisq1D_toy_nh->GetBinCenter(ibin+1),deltachisqval);
        
         chisqval = hchisq1D_toy_ih->GetBinContent(ibin+1);
         hchisq1D_toy_ih_bin[ibin]->Fill(chisqval);
         deltachisqval = hdeltachisq1D_toy_ih->GetBinContent(ibin+1);
         hdeltachisq1D_toy_ih_bin[ibin]->Fill(deltachisqval);
         vec_deltachisq_bin_ih[ibin].push_back(deltachisqval);
        
         hdcp_Vs_deltachisq_ih->Fill(hdeltachisq1D_toy_ih->GetBinCenter(ibin+1),deltachisqval);
        
         vec_like_penalty[ibin].push_back(1.); // no posterior predictive
      }
      
      // Fill histograms of # of events
      double nevents_rhc1re = hErecoToy_RHC1Re->Integral();
      double nevents_1re = hErecoToy_1Re->Integral();
      double nevents_rhc1rmu = hErecoToy_RHC1Rmu->Integral();
      double nevents_1rmu = hErecoToy_1Rmu->Integral();
      
      hNevts_1Rmu->Fill(nevents_1rmu);
      hNevts_1Re->Fill(nevents_1re);
      hNevts_RHC1Rmu->Fill(nevents_rhc1rmu);
      hNevts_RHC1Re->Fill(nevents_rhc1re);
      
   } // end loop over events

   // Scale mean chisq to # of toys and get mean deltachisq (global mean nh and ih)
   hchisq1D_nh_mean->Scale(1./Nentries);
   hchisq1D_ih_mean->Scale(1./Nentries);
    
   TH1D *hdeltachisq1D_nh_mean = (TH1D*) hchisq1D_nh_mean->Clone("hdeltachisq1D_nh_mean");
   TH1D *hdeltachisq1D_ih_mean = (TH1D*) hchisq1D_ih_mean->Clone("hdeltachisq1D_ih_mean");
    
   double glob_min_mean = -999;
   double min_mean_nh = hchisq1D_nh_mean->GetBinContent(hchisq1D_nh_mean->GetMinimumBin());
   double min_mean_ih = hchisq1D_ih_mean->GetBinContent(hchisq1D_ih_mean->GetMinimumBin());
   if(min_mean_nh < min_mean_ih)
   {
      glob_min_mean = min_mean_nh;
   }
   else
   {
      glob_min_mean = min_mean_ih;
   }
    
   for(int ibin=0;ibin<NbinsChisq;ibin++)
   {
      double bincont = hdeltachisq1D_nh_mean->GetBinContent(ibin+1);
      hdeltachisq1D_nh_mean->SetBinContent(ibin+1,bincont-glob_min_mean);
      bincont = hdeltachisq1D_ih_mean->GetBinContent(ibin+1);
      hdeltachisq1D_ih_mean->SetBinContent(ibin+1,bincont-glob_min_mean);
   }
    
   // Build the +-1,2 sigma histograms
    
   TH1D *hchisq1D_nh_p1sigma = hchisq1D_toy_nh->Clone("hchisq1D_nh_p1sigma");
   hchisq1D_nh_p1sigma->Reset();
   TH1D *hchisq1D_ih_p1sigma = hchisq1D_toy_ih->Clone("hchisq1D_ih_p1sigma");
   hchisq1D_ih_p1sigma->Reset();
    
   TH1D *hchisq1D_nh_p2sigma = hchisq1D_toy_nh->Clone("hchisq1D_nh_p1sigma");
   hchisq1D_nh_p2sigma->Reset();
   TH1D *hchisq1D_ih_p2sigma = hchisq1D_toy_ih->Clone("hchisq1D_ih_p1sigma");
   hchisq1D_ih_p2sigma->Reset();
    
   TH1D *hdeltachisq1D_nh_p1sigma = hchisq1D_toy_nh->Clone("hdeltachisq1D_nh_p1sigma");
   hdeltachisq1D_nh_p1sigma->Reset();
   TH1D *hdeltachisq1D_ih_p1sigma = hchisq1D_toy_ih->Clone("hdeltachisq1D_ih_p1sigma");
   hdeltachisq1D_ih_p1sigma->Reset();
    
   TH1D *hdeltachisq1D_nh_p2sigma = hchisq1D_toy_nh->Clone("hdeltachisq1D_nh_p1sigma");
   hdeltachisq1D_nh_p2sigma->Reset();
   TH1D *hdeltachisq1D_ih_p2sigma = hchisq1D_toy_ih->Clone("hdeltachisq1D_ih_p1sigma");
   hdeltachisq1D_ih_p2sigma->Reset();

   // -1,2 sigma histograms
   TH1D *hdeltachisq1D_nh_m1sigma = hchisq1D_toy_nh->Clone("hdeltachisq1D_nh_m1sigma");
   hdeltachisq1D_nh_m1sigma->Reset();
   TH1D *hdeltachisq1D_ih_m1sigma = hchisq1D_toy_ih->Clone("hdeltachisq1D_ih_m1sigma");
   hdeltachisq1D_ih_m1sigma->Reset();
   
   TH1D *hdeltachisq1D_nh_m2sigma = hchisq1D_toy_nh->Clone("hdeltachisq1D_nh_m2sigma");
   hdeltachisq1D_nh_m2sigma->Reset();
   TH1D *hdeltachisq1D_ih_m2sigma = hchisq1D_toy_ih->Clone("hdeltachisq1D_ih_m2sigma");
   hdeltachisq1D_ih_m2sigma->Reset();

   // Median histograms
   TH1D *hdeltachisq1D_nh_median = hchisq1D_toy_nh->Clone("hdeltachisq1D_nh_median");
   hdeltachisq1D_nh_median->Reset();

   TH1D *hdeltachisq1D_ih_median = hchisq1D_toy_ih->Clone("hdeltachisq1D_ih_median");
   hdeltachisq1D_ih_median->Reset();

   // Build p-value histograms
   TH1D *hPvalueChisq_nh = hchisq1D_toy_nh->Clone("hPvalueChisq_nh");
   hPvalueChisq_nh->Reset();
   TH1D *hPvalueChisq_ih = hchisq1D_toy_ih->Clone("hPvalueChisq_ih");
   hPvalueChisq_ih->Reset();
    
   TH1D *hZscoreChisq_nh = hchisq1D_toy_nh->Clone("hZscoreChisq_nh");
   hZscoreChisq_nh->Reset();
   TH1D *hZscoreChisq_ih = hchisq1D_toy_ih->Clone("hZscoreChisq_ih");
   hZscoreChisq_ih->Reset();
    
   TH1D *hProbChisq_90CL_nh = hchisq1D_toy_nh->Clone("hProbChisq_90CL_nh");
   hProbChisq_90CL_nh->Reset();
   TH1D *hProbChisq_90CL_ih = hchisq1D_toy_nh->Clone("hProbChisq_90CL_ih");
   hProbChisq_90CL_ih->Reset();
    
   TH1D *hProbChisq_2sigma_nh = hchisq1D_toy_nh->Clone("hProbChisq_2sigma_nh");
   hProbChisq_2sigma_nh->Reset();
   TH1D *hProbChisq_2sigma_ih = hchisq1D_toy_nh->Clone("hProbChisq_2sigma_ih");
   hProbChisq_2sigma_ih->Reset();

   TH1D *hProbChisq_3sigma_nh = hchisq1D_toy_nh->Clone("hProbChisq_3sigma_nh");
   hProbChisq_3sigma_nh->Reset();
   TH1D *hProbChisq_3sigma_ih = hchisq1D_toy_nh->Clone("hProbChisq_3sigma_ih");
   hProbChisq_3sigma_ih->Reset();

   double maxy_nh = 0;
   double maxy_ih = 0;
   for(int ibin=0;ibin<NbinsChisq;ibin++)
   {
      // Get +1,2sigma critical values of data wrt sensitivity on 1 true point of dcp
      double p1sigma = CalcDLcrit_posterior(vec_deltachisq_bin_nh[ibin],vec_like_penalty[ibin],CLfrac_1sigma*100.,0.);
      double p2sigma = CalcDLcrit_posterior(vec_deltachisq_bin_nh[ibin],vec_like_penalty[ibin],CLfrac_2sigma*100.,0.);
      if(p2sigma > maxy_nh)
      {
         maxy_nh = TMath::Ceil(p2sigma);
      }
      hdeltachisq1D_nh_p1sigma->SetBinContent(ibin+1,p1sigma);
      hdeltachisq1D_nh_p2sigma->SetBinContent(ibin+1,p2sigma);

      p1sigma = CalcDLcrit_posterior(vec_deltachisq_bin_ih[ibin],vec_like_penalty[ibin],CLfrac_1sigma*100.,0.);
      p2sigma = CalcDLcrit_posterior(vec_deltachisq_bin_ih[ibin],vec_like_penalty[ibin],CLfrac_2sigma*100.,0.);
      if(p2sigma > maxy_ih)
      {
         maxy_ih = TMath::Ceil(p2sigma);
      }
      hdeltachisq1D_ih_p1sigma->SetBinContent(ibin+1,p1sigma);
      hdeltachisq1D_ih_p2sigma->SetBinContent(ibin+1,p2sigma);

      // P-value of data wrt sensitivity on 1 true point of dcp
      TH1D *h_nh = (TH1D*) hdeltachisq1D_GlobMin[0][0]; // 0 --> NH
      double chisq_data = h_nh->GetBinContent(ibin+1);
      double pvalue = GetPvalue_rightside(vec_deltachisq_bin_nh[ibin],chisq_data);
      hPvalueChisq_nh->SetBinContent(ibin+1,pvalue);
      TF1 fnh("fnh","ROOT::Math::gaussian_quantile(1.-(x/2.),1.) ",0.,1.);
      double zscore = fnh.Eval(pvalue);
      hZscoreChisq_nh->SetBinContent(ibin+1,zscore);

      TH1D *h_ih = (TH1D*) hdeltachisq1D_GlobMin[0][1]; // 1 --> IH
      chisq_data = h_ih->GetBinContent(ibin+1);
      pvalue = GetPvalue_rightside(vec_deltachisq_bin_ih[ibin],chisq_data);
      hPvalueChisq_ih->SetBinContent(ibin+1,pvalue);
      TF1 fih("fih","ROOT::Math::gaussian_quantile(1.-(x/2.),1.) ",0.,1.);
      zscore = fih.Eval(pvalue);
      hZscoreChisq_ih->SetBinContent(ibin+1,zscore);
      
      // Probability to exclude dcp at 90% or 2sigma CL from sensitivity
      double critval_pt = gFC90CL_nh->Eval(hProbChisq_90CL_nh->GetBinCenter(ibin+1));
      double prob = GetPvalue_rightside(vec_deltachisq_bin_nh[ibin],critval_pt);
      hProbChisq_90CL_nh->SetBinContent(ibin+1,prob);
      critval_pt = gFC2sigma_nh->Eval(hProbChisq_2sigma_nh->GetBinCenter(ibin+1));
      prob = GetPvalue_rightside(vec_deltachisq_bin_nh[ibin],critval_pt);
      hProbChisq_2sigma_nh->SetBinContent(ibin+1,prob);
      critval_pt = gFC3sigma_nh->Eval(hProbChisq_3sigma_nh->GetBinCenter(ibin+1));
      prob = GetPvalue_rightside(vec_deltachisq_bin_nh[ibin],critval_pt);
      hProbChisq_3sigma_nh->SetBinContent(ibin+1,prob);
     
      double critval_pt = gFC90CL_ih->Eval(hProbChisq_90CL_ih->GetBinCenter(ibin+1));
      double prob = GetPvalue_rightside(vec_deltachisq_bin_ih[ibin],critval_pt);
      hProbChisq_90CL_ih->SetBinContent(ibin+1,prob);
      critval_pt = gFC2sigma_ih->Eval(hProbChisq_2sigma_ih->GetBinCenter(ibin+1));
      prob = GetPvalue_rightside(vec_deltachisq_bin_ih[ibin],critval_pt);
      hProbChisq_2sigma_ih->SetBinContent(ibin+1,prob);
      critval_pt = gFC3sigma_ih->Eval(hProbChisq_3sigma_ih->GetBinCenter(ibin+1));
      prob = GetPvalue_rightside(vec_deltachisq_bin_ih[ibin],critval_pt);
      hProbChisq_3sigma_ih->SetBinContent(ibin+1,prob);
   }

   // Print
   int bin_mix23_nova = hProbChisq_2sigma_ih->FindBin(exclusion_value);
   double prob_excl_dcp0_90CL_nh = hProbChisq_90CL_nh->GetBinContent(bin_mix23_nova);
   double prob_excl_dcp0_90CL_ih = hProbChisq_90CL_ih->GetBinContent(bin_mix23_nova);
   double prob_excl_dcp0_2sigma_nh = hProbChisq_2sigma_nh->GetBinContent(bin_mix23_nova);
   double prob_excl_dcp0_2sigma_ih = hProbChisq_2sigma_ih->GetBinContent(bin_mix23_nova);
   double prob_excl_dcp0_3sigma_nh = hProbChisq_3sigma_nh->GetBinContent(bin_mix23_nova);
   double prob_excl_dcp0_3sigma_ih = hProbChisq_3sigma_ih->GetBinContent(bin_mix23_nova);
  
   cout << "Normal Hierarchy:" << endl;
   cout << " - Prob. to exclude T2K bestfit at 90%CL: " << prob_excl_dcp0_90CL_nh << endl;
   cout << " - Prob. to exclude T2K bestfit at 2sigma: " << prob_excl_dcp0_2sigma_nh << endl;
   cout << " - Prob. to exclude T2K bestfit at 3sigma: " << prob_excl_dcp0_3sigma_nh << endl;
   cout << "Inverted Hierarchy:" << endl;
   cout << " - Prob. to exclude T2K bestfit at 90%CL: " << prob_excl_dcp0_90CL_ih << endl;
   cout << " - Prob. to exclude T2K bestfit at 2sigma: " << prob_excl_dcp0_2sigma_ih << endl;
   cout << " - Prob. to exclude T2K bestfit at 3sigma: " << prob_excl_dcp0_3sigma_ih << endl;
}

double CalcBestFitErrorsMINOS(TH1D *histo, double CritVal, double val_bf){
  
  cout << endl;
  cout << "Compute MINOS errors" << endl;
  cout << endl;
  
  //double error = conv::kBadResult;
  const double CritChisq_1D = CritVal;
  
  // Find 68.27% CL contour (1sigma)
  
  double errneg = conv::kBadResult;
  double errpos = conv::kBadResult;
  
  int bin_bf = conv::kBadResult;        // best-fit bin
  double chisq_min = conv::kBadResult;  // best-fit chisq
  
  if(val_bf==conv::kBadResult){
    chisq_min = histo->GetBinContent(histo->GetMinimumBin());
    bin_bf = histo->GetMinimumBin();
    val_bf = histo->GetBinCenter(bin_bf);
    
    //cout << endl;
    cout << "Find the global minimum:" << endl;
    cout << "best-fit value = " << val_bf << endl;
    cout << "best-fit chisq = " << chisq_min << endl;
    //cout << endl;
  }
  else{
    bin_bf = histo->FindBin(val_bf);             // best-fit bin
    chisq_min = histo->GetBinContent(bin_bf); // best-fit chisq
    
    //cout << endl;
    cout << "Find the global minimum:" << endl;
    cout << "best-fit value = " << val_bf << endl;
    cout << "best-fit chisq = " << chisq_min << endl;
    //cout << endl;
  }
  
  
  // cout << endl;
  // cout << "val_bf = " << val_bf << endl;
  // cout << "bin_bf = " << bin_bf << endl;
  // cout << "chisq_min = " << chisq_min << endl;
  // cout << endl;
  
  // Bins before (b) and after (a) the intersection with critical value
  int bin_negsigma_a = conv::kBadResult;
  int bin_negsigma_b = conv::kBadResult;
  int bin_possigma_a = conv::kBadResult;
  int bin_possigma_b = conv::kBadResult;
  
  // Find param value at DeltaChisq=1 (negative error)
  
  for(int ibin=bin_bf;ibin>=0;ibin--){ // histogram convention
    
    double bincont = histo->GetBinContent(ibin);
    // if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Bin content is negative!" << endl;
    //   cerr << "bincont = " << bincont << endl;
    //   cerr << endl;
    // }
    
    // cout << "ibin = " << ibin
    // 	   << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
    // 	   << " - bincont = " << bincont << endl;
    
    if(bincont>CritChisq_1D){
      // bins before (b) and after (a) DeltaChisq=1 (following x direction)
      bin_negsigma_b = ibin;
      bin_negsigma_a = bin_negsigma_b+1;
      break;
    }
    // else if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Negative DeltaChisq 1D:" << endl;
    //   cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
    //   cerr << endl;
    //   //exit(1);
    //   return conv::kBadResult;
    // }
    
  }
  
  
  // Positive error
  
  // cout << endl;
  // cout << "Positive error" << endl;
  // cout << endl;
  
  for(int ibin=bin_bf;ibin<=histo->GetNbinsX();ibin++){ // histogram convention
    double bincont = histo->GetBinContent(ibin);
    // if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Negative DeltaChisq 1D:" << endl;
    //   cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
    //   cerr << endl;
    //   //exit(1);
    //   return conv::kBadResult;
    // }
    
    // cout << "ibin = " << ibin
    // 	   << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
    //  	   << " - bincont = " << bincont << endl;
    
    if(bincont>CritChisq_1D){
      // bins before (b) and after (a) DeltaChisq=1
      bin_possigma_a = ibin;
      bin_possigma_b = bin_possigma_a-1;
      break;
    }
  }
  //cout << "bin_possigma_b = " << bin_possigma_b << endl;
  //cout << "bin_possigma_a = " << bin_possigma_a << endl;
  
  
  // If no bin is found it means the 1sigma bin is outside the physical region
  // and bin variables are still set to default value.
  // If so set left extreme to 1 (first bin of the histogram) and right extreme
  // to Nbins (last bin of the histogram)
  
  if(bin_negsigma_a == conv::kBadResult){
    bin_negsigma_a = 1;
    bin_negsigma_b = 1;
  }
  if(bin_possigma_b == conv::kBadResult){
    bin_possigma_a = histo->GetNbinsX();
    bin_possigma_b = histo->GetNbinsX();
  }
  
  // Check if the interval bins are correctly found
  if(bin_negsigma_a > bin_possigma_b){
    cerr << endl;
    cerr << "FitMarginalization::CalcBestFitErrorsMINOS" << endl;
    cerr << "Left extreme bin is higher than right extreme bin" << endl;
    cerr << "The 1sigma error interval is messed up:" << endl;
    cerr << "bin_possigma_b = " << bin_possigma_b << endl;
    cerr << "bin_negsigma_a = " << bin_negsigma_a << endl;
    cerr << endl;
    //exit(1);
    return conv::kBadResult;
  }
  
  double X1 = conv::kBadResult;
  double Y1 = conv::kBadResult;
  double X2 = conv::kBadResult;
  double Y2 = conv::kBadResult;
  double Xout = conv::kBadResult;
  
  // Extrapolation (left extreme). Do it only if 1sigma bin is in the physical region
  if(bin_negsigma_a>1){
    X1 = histo->GetXaxis()->GetBinCenter(bin_negsigma_b);
    Y1 = histo->GetBinContent(bin_negsigma_b);
    X2 = histo->GetXaxis()->GetBinCenter(bin_negsigma_a);
    Y2 = histo->GetBinContent(bin_negsigma_a);
    Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
  }
  else{
    Xout = histo->GetBinCenter(bin_negsigma_a);
  }
  
  errneg = fabs(val_bf - Xout);
  //cout << "fparam_bf_errneg[0] = " << fparam_bf_errneg[0] << endl;
  
  // Extrapolation (right extreme). Do it only if 1sigma bin is in the physical region
  if(bin_possigma_b < histo->GetNbinsX()){
    X1 = histo->GetXaxis()->GetBinCenter(bin_possigma_b);
    Y1 = histo->GetBinContent(bin_possigma_b);
    X2 = histo->GetXaxis()->GetBinCenter(bin_possigma_a);
    Y2 = histo->GetBinContent(bin_possigma_a);
    Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
  }
  else{
    Xout = histo->GetBinCenter(bin_possigma_b);
  }
  errpos[0] = fabs(val_bf - Xout);
  //cout << "fparam_bf_errpos[0] = " << fparam_bf_errpos[0] << endl;
  
  cout << endl;
  cout << "Best-fit = " << val_bf << " + " << errpos << " - " << errneg << endl;
  cout << endl;
  cout << "1 sigma range:"<< endl;
  cout << "[ " << val_bf - errneg << " ; " << val_bf + errpos << " ] " << endl;
  cout << endl;
}


double CalcBestFitErrorsMINOS_FC(TH1D *histo, double CritVal, double val_bf,double &lim_x_min, double &lim_x_max){
  
  cout << endl;
  cout << "Feldman-Cousins confidence intervals" << endl;
  cout << endl;
  
  //double error = conv::kBadResult;
  const double CritChisq_1D = CritVal;
  
  // Find 68.27% CL contour (1sigma)
  
  double errneg = conv::kBadResult;
  double errpos = conv::kBadResult;
  
  int bin_bf = conv::kBadResult;        // best-fit bin
  double chisq_min = conv::kBadResult;  // best-fit chisq
  
  if(val_bf==conv::kBadResult){
    chisq_min = histo->GetBinContent(histo->GetMinimumBin());
    bin_bf = histo->GetMinimumBin();
    val_bf = histo->GetBinCenter(bin_bf);
    
    //cout << endl;
    cout << "Find the global minimum:" << endl;
    cout << "best-fit value = " << val_bf << endl;
    cout << "best-fit chisq = " << chisq_min << endl;
    //cout << endl;
  }
  else{
    bin_bf = histo->FindBin(val_bf);             // best-fit bin
    chisq_min = histo->GetBinContent(bin_bf); // best-fit chisq
    
    //cout << endl;
    cout << "Find the global minimum:" << endl;
    cout << "best-fit value = " << val_bf << endl;
    cout << "best-fit chisq = " << chisq_min << endl;
    //cout << endl;
  }
  
  
  // cout << endl;
  // cout << "val_bf = " << val_bf << endl;
  // cout << "bin_bf = " << bin_bf << endl;
  // cout << "chisq_min = " << chisq_min << endl;
  // cout << endl;
  
  // Bins before (b) and after (a) the intersection with critical value
  int bin_negsigma_a = conv::kBadResult;
  int bin_negsigma_b = conv::kBadResult;
  int bin_possigma_a = conv::kBadResult;
  int bin_possigma_b = conv::kBadResult;
  
  // Find param value at DeltaChisq=1 (negative error)
  
  for(int ibin=bin_bf;ibin>=0;ibin--){ // histogram convention
    
    double bincont = histo->GetBinContent(ibin);
    // if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Bin content is negative!" << endl;
    //   cerr << "bincont = " << bincont << endl;
    //   cerr << endl;
    // }
    
    //
    // Useful for checks
    //
    // cout << "ibin = " << ibin
    // 	   << " - X = " << histo->GetXaxis()->GetBinCenter(ibin)
    // 	   << " - bincont = " << bincont << endl;
    
    
    if(bincont>CritChisq_1D){
      // bins before (b) and after (a) DeltaChisq=1 (following x direction)
      bin_negsigma_b = ibin;
      bin_negsigma_a = bin_negsigma_b+1;
      break;
    }
    // else if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Negative DeltaChisq 1D:" << endl;
    //   cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
    //   cerr << endl;
    //   //exit(1);
    //   return conv::kBadResult;
    // }
    
  }
  
  
  //
  // Useful for checks
  //
  // double x1 = histo->GetBinCenter(bin_negsigma_b);
  // double y1 = histo->GetBinContent(bin_negsigma_b);
  // double x2 = histo->GetBinCenter(bin_negsigma_a);
  // double y2 = histo->GetBinContent(bin_negsigma_a);
  // cout << "x1 = " << x1 << ", y1 = " << y1 << endl;
  // cout << "x2 = " << x2 << ", y2 = " << y2 << endl;
  // cout << " - interpolation = " << GetInterpolationX(x1,y1,x2,y2,0.) << endl;
  
  
  
  
  // Positive error
  
  // cout << endl;
  // cout << "Positive error" << endl;
  // cout << endl;
  
  for(int ibin=bin_bf;ibin<=histo->GetNbinsX();ibin++){ // histogram convention
    double bincont = histo->GetBinContent(ibin);
    // if(bincont<0.){
    //   cerr << endl;
    //   cerr << "CalcBestFitErrorsMINOS" << endl;
    //   cerr << "Negative DeltaChisq 1D:" << endl;
    //   cerr << "DeltaChisq_1D = " << histo->GetBinContent(ibin) << endl;
    //   cerr << endl;
    //   //exit(1);
    //   return conv::kBadResult;
    // }
    
    // cout << "ibin = " << ibin
    // 	   << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
    //  	   << " - bincont = " << bincont << endl;
    
    if(bincont>CritChisq_1D){
      // bins before (b) and after (a) DeltaChisq=1
      bin_possigma_a = ibin;
      bin_possigma_b = bin_possigma_a-1;
      break;
    }
  }
  //cout << "bin_possigma_b = " << bin_possigma_b << endl;
  //cout << "bin_possigma_a = " << bin_possigma_a << endl;
  
  
  // If no bin is found it means the 1sigma bin is outside the physical region
  // and bin variables are still set to default value.
  // If so set left extreme to 1 (first bin of the histogram) and right extreme
  // to Nbins (last bin of the histogram)
  
  if(bin_negsigma_a == conv::kBadResult){
    bin_negsigma_a = 1;
    bin_negsigma_b = 1;
  }
  if(bin_possigma_b == conv::kBadResult){
    bin_possigma_a = histo->GetNbinsX();
    bin_possigma_b = histo->GetNbinsX();
  }
  
  // Check if the interval bins are correctly found
  if(bin_negsigma_a > bin_possigma_b){
    cerr << endl;
    cerr << "FitMarginalization::CalcBestFitErrorsMINOS" << endl;
    cerr << "Left extreme bin is higher than right extreme bin" << endl;
    cerr << "The 1sigma error interval is messed up:" << endl;
    cerr << "bin_possigma_b = " << bin_possigma_b << endl;
    cerr << "bin_negsigma_a = " << bin_negsigma_a << endl;
    cerr << endl;
    //exit(1);
    return conv::kBadResult;
  }
  
  double X1 = conv::kBadResult;
  double Y1 = conv::kBadResult;
  double X2 = conv::kBadResult;
  double Y2 = conv::kBadResult;
  double Xout = conv::kBadResult;
  
  // Extrapolation (left extreme). Do it only if 1sigma bin is in the physical region
  if(bin_negsigma_a>1){
    X1 = histo->GetXaxis()->GetBinCenter(bin_negsigma_b);
    Y1 = histo->GetBinContent(bin_negsigma_b);
    X2 = histo->GetXaxis()->GetBinCenter(bin_negsigma_a);
    Y2 = histo->GetBinContent(bin_negsigma_a);
    Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
  }
  else{
    Xout = histo->GetBinCenter(bin_negsigma_a);
  }
  
  errneg = fabs(val_bf - Xout);
  //cout << "fparam_bf_errneg[0] = " << fparam_bf_errneg[0] << endl;
  
  // Extrapolation (right extreme). Do it only if 1sigma bin is in the physical region
  if(bin_possigma_b < histo->GetNbinsX()){
    X1 = histo->GetXaxis()->GetBinCenter(bin_possigma_b);
    Y1 = histo->GetBinContent(bin_possigma_b);
    X2 = histo->GetXaxis()->GetBinCenter(bin_possigma_a);
    Y2 = histo->GetBinContent(bin_possigma_a);
    Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
  }
  else{
    Xout = histo->GetBinCenter(bin_possigma_b);
  }
  errpos[0] = fabs(val_bf - Xout);
  //cout << "fparam_bf_errpos[0] = " << fparam_bf_errpos[0] << endl;
  
  cout << endl;
  cout << "Best-fit = " << val_bf << " + " << errpos << " - " << errneg << endl;
  cout << endl;
  cout << "The confidence interval is:"<< endl;
  cout << "[ " << val_bf - errneg << " ; " << val_bf + errpos << " ] " << endl;
  cout << endl;
  
  lim_x_min = val_bf - errneg;
  lim_x_max = val_bf + errpos;
  
}



double GetInterpolationX(double X1,double Y1,double X2, double Y2,double Yin){
  
  double slope = (Y2-Y1)/(X2-X1);
  double Xout = (Yin-Y1)/slope + X1;
  return Xout;
  
  // cout << "Yin = " << Yin << endl;
  // cout << "X1 = " << X1 << " - X2 = " << X2 << endl;
  // cout << "Y1 = " << Y1 << " - Y2 = " << Y2 << endl;
  // cout << "Xout = " << Xout << endl;
  // cout << endl;
  // double Xin = Xout;
  // Yin = (Y2-Y1)/(X2-X1)*(Xin-X1) + Y1;
  // cout << "Yin = " << Yin << endl;
  
  // TGraph *gr = new TGraph();
  // gr->SetPoint(0,X1,Y1);
  // gr->SetPoint(1,X2,Y2);
  // cout << "Eval = " << gr->Eval(Xin) << endl;
  // return gr->Eval(Xin);
}


double CalcDLcrit_posterior(vector<double> dchisq, vector<double> weight, double conflev, double sigmacl){
  
  // Get error on coverage: sigmacl = sqrt(N*p*(1-p))/N
  
  int nelem = dchisq.size();
  
  double *new_dchisq = new double[nelem];
  double *new_weight = new double[nelem];
  int *index = new int[nelem];
  
  // copy vector into array
  for(int i=0;i<nelem;i++){
    new_dchisq[i] = dchisq[i];
    new_weight[i] = weight[i];
    index[i] = i;
  }
  
  sort(dchisq.begin(), dchisq.end());
  
  // store the sorted indices --> use for both dchisq and weights
  //index = TMath::Sort(nelem,new_dchisq,index,0);
  TMath::Sort(nelem,new_dchisq,index,0);
  
  double tot_weight = 0.;
  for(int i=0;i<nelem;i++) tot_weight += new_weight[index[i]];
  
  //cout << "conflev = " << conflev << endl;
  
  double sum = 0.;
  double critval = -999.;
  for(int i=0;i<nelem;i++){
    sum += new_weight[index[i]];
    
    //cout << std::setprecision(16);
    //cout << sum << " ";
    
    double epsilon = 0.000001; // for the double precision in the break condition

    if( sum >= ( (0.01 * conflev + sigmacl)*tot_weight - epsilon ) ){
      
      // cout << endl;
      // cout << "sum = " << sum << endl;
      // cout << "new_weight[index[i]] = " << new_weight[index[i]] << endl;
      // cout << "before: sum = " << (sum - new_weight[index[i]]) << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      // cout << "current: sum = " << sum << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      // cout << "after: sum = " << (sum + new_weight[index[i]]) << " - CL = " << (0.01 * conflev + sigmacl)*tot_weight << endl;
      
      critval = new_dchisq[index[i]];
      break;
    }
  }
  
  return critval;
}

std::pair<double, double> CalcDLcrit_binned(vector<double> dchisq, double conflev, double& median)
{
  int nelem = dchisq.size();
  double cut_off = conflev * nelem;
  double min_dchisq = *std::max_element(dchisq.begin(), dchisq.end());
  double max_dchisq = *std::min_element(dchisq.begin(), dchisq.end());
  const int nbins = 100;
  TH1D h_dchisq = TH1D("dchisq", "dchisq", nbins, min_dchisq - 1e-6, max_dchisq + 1e-6);
  for(int i = 0; i < nelem; ++i)
  {
     h_dchisq.Fill(dchisq[i]);
  }
  int* toy_count = new int[nbins];
  for(int i = 0; i < nbins; ++i)
  {
     toy_count[i] = h_dchisq.GetBinContent(i + 1);
  }
  int* index = new int[nbins];
  
  // store the sorted indices in descending order
  TMath::Sort(nbins, toy_count, index, true);

  int interval_size = 0;

  for(int i = 0; i < nbins; ++i)
  {
     interval_size += toy_count[index[i]];
     double low_dchisq = h_dchisq.GetXaxis()->GetBinLowEdge(index[i] + 1);
     double upp_dchisq = h_dchisq.GetXaxis()->GetBinUpEdge(index[i] + 1);
     if(low_dchisq < min_dchisq)
     {
        min_dchisq = low_dchisq;
     }
     if(upp_dchisq > max_dchisq)
     {
        max_dchisq = upp_dchisq;
     }
     if(interval_size >= cut_off)
     {
        break;
     }
  }

  std::sort(dchisq.begin(), dchisq.end());
  int idx = nelem / 2;
  if((nelem % 2) == 0)
  {
     median = (dchisq[idx - 1] + dchisq[idx]) / 2.0;
  }
  else
  {
     median = dchisq[idx];
  }

  return std::make_pair(min_dchisq, max_dchisq);
}


double CalcDLcrit(vector<double> dchisq, double conflev, double sigmacl){
  
  // Get error on critical value: sigmacl = sqrt(N*p*(1-p))
  
  sort(dchisq.begin(), dchisq.end()); // sort the dchisq vector
  vector<double>::iterator dchisqIter;
  dchisqIter = dchisq.begin() + (int)(0.01*conflev*dchisq.size() + sigmacl);
  double critval = (*dchisqIter);
  return critval;
}

double CalcDLcrit(TH1D *h_delta_true,double FCfrac_CL){
  
  //Feldman-Cousins: Calculate chi2 value at 90%CL************************
  double ntoys_total = 0;
  double ntoys = 0;
  double ntoys_frac = 0;
  int chi2_bin = 0;
  double chi2_val = 0;
  double dlval = 0;
  
  ntoys_total = h_delta_true->Integral();
  
  for(int ibin=0;ibin<h_delta_true->GetNbinsX();ibin++){
    ntoys += h_delta_true->GetBinContent(ibin+1);
    ntoys_frac = ntoys/ntoys_total;
    
    if(ntoys_frac>=FCfrac_CL) {
      chi2_val = h_delta_true->GetBinCenter(ibin+1);
      break;
    }
  }
  return chi2_val;
}


double GetPvalue_rightside(vector<double> vec_toys,double data){
  
  double ntot_toys = vec_toys.size();
  
  double nfrac=0;
  for (std::vector<double>::iterator it=vec_toys.begin(); it!=vec_toys.end(); ++it){
    if(*it >= data) nfrac += 1.; // If toy=data include the entry in the pvalue integral range
  }
  double pvalue =  nfrac / (double)ntot_toys;
  
  return pvalue;
}

double GetJointPvalue_rightside(vector<double> vec_0_toys, vector<double> vec_pi_toys,
      double critical_0, double critical_pi)
{
   int ntot_toys = vec_0_toys.size();
   double nfrac = 0;
   for (int i = 0; i < ntot_toys; ++i)
   {
      double dchisq0 = vec_0_toys[i];
      double dchisqpi = vec_pi_toys[i];
      if(dchisq0 >= critical_0 && dchisqpi >= critical_pi)
      {
         nfrac += 1.;
      }
   }
   double pvalue =  nfrac / ntot_toys;
   return pvalue;
}


double GetPvalue_posterior_rightside(vector<double> vec_toys,vector<double> weight,double data){
  
  double ntot_toys = vec_toys.size();
  
  double tot = 0.;
  double nfrac=0;
  for(int i=0;i<vec_toys.size();i++){
    tot += weight[i];
    double cont = vec_toys[i];
    if(cont >= data) nfrac += weight[i];
  }
  double pvalue =  nfrac / tot;
  
  return pvalue;
}

struct fit_info {
  int bf_hierarchy;
  double dchisq;
  double bf_dcp;
  
//  bool operator<( const fit_info& val ) const {
//    return dchisq < val.dchisq;
//  }

  fit_info(int _bf_hierarchy, double _dchisq, double _bf_dcp) : bf_hierarchy(_bf_hierarchy), dchisq(_dchisq), bf_dcp(_bf_dcp) {}

};

//struct sort_by_dchisq
//{
//  inline bool operator() (const fit_info& struct1, const fit_info& struct2)
//  {
//    return (struct1.dchisq < struct2.dchisq);
//  }
//};

int format_hist(TH1D* _histogram)
{
  _histogram->GetYaxis()->SetTitleOffset(0.85);
  _histogram->GetXaxis()->SetTitleOffset(1);
  _histogram->GetXaxis()->SetLabelSize(0.05);
  _histogram->GetYaxis()->SetLabelSize(0.05);
  _histogram->GetXaxis()->SetLabelFont(132);
  _histogram->GetYaxis()->SetLabelFont(132);
  _histogram->GetXaxis()->SetTitleFont(132);
  _histogram->GetYaxis()->SetTitleFont(132);
  _histogram->GetXaxis()->SetTitleSize(0.06);
  _histogram->GetYaxis()->SetTitleSize(0.06);
  _histogram->GetXaxis()->CenterTitle();
  _histogram->GetYaxis()->CenterTitle();
  _histogram->GetZaxis()->CenterTitle();
}


void plot_slice(vector<double> dchisq, double data_chisq, double dcp_slice, const char* filename)
{ 
   int nelem = dchisq.size();
   double* new_dchisq = new double[nelem];
  
   sort(dchisq.begin(), dchisq.end());
  
   double dchisq_max = dchisq[nelem - 1];
  
   TH1D dchisq_dist("hdchisq_slice", Form("Slice through #delta_{CP} = %.2f", dcp_slice), 100, 0.0, dchisq_max);
   dchisq_dist.SetLineColor(kBlack);
   int cutoff = -1;
   for(int i = 0; i < nelem;i++)
   {
      dchisq_dist.Fill(dchisq[i]);
      if(dchisq[i] < data_chisq)
      {
         cutoff = i;
      }
   }
   double pvalue = 1.0 - (cutoff + 1.0) / nelem;
   dchisq_dist.Scale(1.0 / dchisq_dist.Integral());
   dchisq_dist.GetXaxis()->SetTitle("-2#Deltaln(L)");
   dchisq_dist.GetYaxis()->SetTitle("Fraction of toys");

   TCanvas c_slice("c_slice","c_slice");
   c_slice.SetTopMargin(0.1);
   dchisq_dist.Draw("hist");
   TLine ldata(data_chisq, 0, data_chisq, dchisq_dist.GetMaximum());
   ldata.SetLineColor(kRed);
   ldata.SetLineWidth(2);
   ldata.Draw("same");
   TLatex title;
   title.SetTextSize(0.04);
   double x = dchisq_dist.GetBinCenter((int)(dchisq_dist.GetNbinsX() * 0.5));
   double y = dchisq_dist.GetMaximum() * 1.1;
   title.DrawLatex(x, y, dchisq_dist.GetTitle());
   TLatex text;
   text.SetTextSize(0.04);
   x = dchisq_dist.GetBinCenter((int)(dchisq_dist.GetNbinsX() * 0.75));
   y = dchisq_dist.GetMaximum() * 0.75;
   text.DrawLatex(x, y, Form("p-value = %.3f", pvalue));
   c_slice.Print(filename);
}

