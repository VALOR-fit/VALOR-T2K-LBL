//
// ------------------------------------------------------------------------------------------------
// Plot 1-ring mu E_{reco} spectrum.
// The spectra can be rebinned and the different modes can be regrouped.
// Data and alternative spectra can be superimposed.
// Works with the output of the macros/spectra/test_sk_ereco_spectrum_*.C set of macros.
//
// C.Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void plot_1ringe_spectrum_rebin_regroup
(
  /* ROOT file with pdf histograms (total + contributions from various modes) */
  const char * pdffilename,                 
  const char * pdflabel = "pdf", /* legend label for pdf, "best-fit, Run1+2+3" */

  /* (optional) ROOT file with alternative pdf histograms (total + contributions from various modes) */
  const char * altpdffilename = 0,                 
  const char * altpdflabel    = "",  /* legend label for alternative pdf */

  /* (optional) ROOT file with data to superimpose */
  const char * datafilename  = 0,                 
  const char * datahistname  = "",  /* data histogram name */
  const char * datalabel     = "",  /* data label */

  /* re-bin E_{reco} spectra (fewer variable size bins) */  
  bool rebin = false,    

  /* max E_{reco} in plot  */
  double emax = 1.25,      

  /* max of N_{events} axis in plot  */
  double nmax = 200.0,

  /* show pdf components (shown only for the first pdf, for clarity) */
  bool show_components = true,

  /* use "1 #mu-like ring" as a legend header */
  bool print_header = false
)
{
  // 
  // Alternative binning scheme
  //
  const int kNErecoBins1Rmu = 16;
  const double kErecoBinEdges1Rmu[kNErecoBins1Rmu+1] = {
   0.000,  0.250,  0.400,  0.500,
   0.600,  0.700,  0.800,  0.900,
   1.000,  1.100,  1.250,  1.500,  2.000, 3.000,  4.000, 6.000, 20.000
  };

  //
  // Get pdfs (total & components) from input file
  //

  TFile * pdffile  = new TFile(pdffilename, "read");
  pdffile->cd();
  TH1D * sk_1Retotal          = (TH1D*) pdffile->Get("sk_1Retotal");
  TH1D * sk_1Re_numuccqe      = (TH1D*) pdffile->Get("sk_1Re_numuccqe");
  TH1D * sk_1Re_numucc1pi     = (TH1D*) pdffile->Get("sk_1Re_numucc1pi");
  TH1D * sk_1Re_numucccoh     = (TH1D*) pdffile->Get("sk_1Re_numucccoh");
  TH1D * sk_1Re_numuccoth     = (TH1D*) pdffile->Get("sk_1Re_numuccoth");
  TH1D * sk_1Re_numunc1pi0    = (TH1D*) pdffile->Get("sk_1Re_numunc1pi0");
  TH1D * sk_1Re_numunccoh    = (TH1D*) pdffile->Get("sk_1Re_numunccoh");
  TH1D * sk_1Re_numuncoth     = (TH1D*) pdffile->Get("sk_1Re_numuncoth");
  TH1D * sk_1Re_oscnumuccqe   = (TH1D*) pdffile->Get("sk_1Re_oscnumuccqe");
  TH1D * sk_1Re_oscnumucc1pi  = (TH1D*) pdffile->Get("sk_1Re_oscnumucc1pi");
  TH1D * sk_1Re_oscnumucccoh  = (TH1D*) pdffile->Get("sk_1Re_oscnumucccoh");
  TH1D * sk_1Re_oscnumuccoth  = (TH1D*) pdffile->Get("sk_1Re_oscnumuccoth");
  TH1D * sk_1Re_numubarccqe   = (TH1D*) pdffile->Get("sk_1Re_numubarccqe");
  TH1D * sk_1Re_numubarcc1pi  = (TH1D*) pdffile->Get("sk_1Re_numubarcc1pi");
  TH1D * sk_1Re_numubarcccoh  = (TH1D*) pdffile->Get("sk_1Re_numubarcccoh");
  TH1D * sk_1Re_numubarccoth  = (TH1D*) pdffile->Get("sk_1Re_numubarccoth");
  TH1D * sk_1Re_numubarnc1pi0 = (TH1D*) pdffile->Get("sk_1Re_numubarnc1pi0");
  TH1D * sk_1Re_numubarnccoh = (TH1D*) pdffile->Get("sk_1Re_numubarnccoh");
  TH1D * sk_1Re_numubarncoth  = (TH1D*) pdffile->Get("sk_1Re_numubarncoth");
  TH1D * sk_1Re_nueccqe       = (TH1D*) pdffile->Get("sk_1Re_nueccqe");
  TH1D * sk_1Re_nuecc1pi      = (TH1D*) pdffile->Get("sk_1Re_nuecc1pi");
  TH1D * sk_1Re_nuecccoh      = (TH1D*) pdffile->Get("sk_1Re_nuecccoh");
  TH1D * sk_1Re_nueccoth      = (TH1D*) pdffile->Get("sk_1Re_nueccoth");
  TH1D * sk_1Re_nuenc1pi0     = (TH1D*) pdffile->Get("sk_1Re_nuenc1pi0");
  TH1D * sk_1Re_nuenccoh     = (TH1D*) pdffile->Get("sk_1Re_nuenccoh");
  TH1D * sk_1Re_nuencoth      = (TH1D*) pdffile->Get("sk_1Re_nuencoth");
  TH1D * sk_1Re_nuebarccqe    = (TH1D*) pdffile->Get("sk_1Re_nuebarccqe");
  TH1D * sk_1Re_nuebarcc1pi   = (TH1D*) pdffile->Get("sk_1Re_nuebarcc1pi");
  TH1D * sk_1Re_nuebarcccoh   = (TH1D*) pdffile->Get("sk_1Re_nuebarcccoh");
  TH1D * sk_1Re_nuebarccoth   = (TH1D*) pdffile->Get("sk_1Re_nuebarccoth");
  TH1D * sk_1Re_nuebarnc1pi0  = (TH1D*) pdffile->Get("sk_1Re_nuebarnc1pi0");
  TH1D * sk_1Re_nuebarnccoh  = (TH1D*) pdffile->Get("sk_1Re_nuebarnccoh");
  TH1D * sk_1Re_nuebarncoth   = (TH1D*) pdffile->Get("sk_1Re_nuebarncoth");
  TH1D * sk_1Re_oscnumubarccqe    = (TH1D*) pdffile->Get("sk_1Re_oscnumubarccqe");
  TH1D * sk_1Re_oscnumubarcc1pi   = (TH1D*) pdffile->Get("sk_1Re_oscnumubarcc1pi");
  TH1D * sk_1Re_oscnumubarcccoh   = (TH1D*) pdffile->Get("sk_1Re_oscnumubarcccoh");
  TH1D * sk_1Re_oscnumubarccoth   = (TH1D*) pdffile->Get("sk_1Re_oscnumubarccoth");
  TH1D * sk_1Re_oscnueccqe    = (TH1D*) pdffile->Get("sk_1Re_oscnueccqe");
  TH1D * sk_1Re_oscnuecc1pi   = (TH1D*) pdffile->Get("sk_1Re_oscnuecc1pi");
  TH1D * sk_1Re_oscnuecccoh   = (TH1D*) pdffile->Get("sk_1Re_oscnuecccoh");
  TH1D * sk_1Re_oscnueccoth   = (TH1D*) pdffile->Get("sk_1Re_oscnueccoth");
  TH1D * sk_1Re_oscnuebarccqe    = (TH1D*) pdffile->Get("sk_1Re_oscnuebarccqe");
  TH1D * sk_1Re_oscnuebarcc1pi   = (TH1D*) pdffile->Get("sk_1Re_oscnuebarcc1pi");
  TH1D * sk_1Re_oscnuebarcccoh   = (TH1D*) pdffile->Get("sk_1Re_oscnuebarcccoh");
  TH1D * sk_1Re_oscnuebarccoth   = (TH1D*) pdffile->Get("sk_1Re_oscnuebarccoth");

  if(!sk_1Retotal) sk_1Retotal = (TH1D*) pdffile->Get("sk_RHC1Retotal");
  if(!sk_1Re_numuccqe) sk_1Re_numuccqe = (TH1D*) pdffile->Get("sk_RHC1Re_numuccqe");
  if(!sk_1Re_numucc1pi) sk_1Re_numucc1pi = (TH1D*) pdffile->Get("sk_RHC1Re_numucc1pi");
  if(!sk_1Re_numucccoh) sk_1Re_numucccoh = (TH1D*) pdffile->Get("sk_RHC1Re_numucccoh");
  if(!sk_1Re_numuccoth) sk_1Re_numuccoth = (TH1D*) pdffile->Get("sk_RHC1Re_numuccoth");
  if(!sk_1Re_numunc1pi0) sk_1Re_numunc1pi0 = (TH1D*) pdffile->Get("sk_RHC1Re_numunc1pi0");
  if(!sk_1Re_numunccoh) sk_1Re_numunccoh = (TH1D*) pdffile->Get("sk_RHC1Re_numunccoh");
  if(!sk_1Re_numuncoth) sk_1Re_numuncoth = (TH1D*) pdffile->Get("sk_RHC1Re_numuncoth");
  if(!sk_1Re_oscnumuccqe) sk_1Re_oscnumuccqe = (TH1D*) pdffile->Get("sk_RHC1Re_oscnumuccqe");
  if(!sk_1Re_oscnumucc1pi) sk_1Re_oscnumucc1pi = (TH1D*) pdffile->Get("sk_RHC1Re_oscnumucc1pi");
  if(!sk_1Re_oscnumucccoh) sk_1Re_oscnumucccoh = (TH1D*) pdffile->Get("sk_RHC1Re_oscnumucccoh");
  if(!sk_1Re_oscnumuccoth) sk_1Re_oscnumuccoth = (TH1D*) pdffile->Get("sk_RHC1Re_oscnumuccoth");
  if(!sk_1Re_numubarccqe) sk_1Re_numubarccqe = (TH1D*) pdffile->Get("sk_RHC1Re_numubarccqe");
  if(!sk_1Re_numubarcc1pi) sk_1Re_numubarcc1pi = (TH1D*) pdffile->Get("sk_RHC1Re_numubarcc1pi");
  if(!sk_1Re_numubarcccoh) sk_1Re_numubarcccoh = (TH1D*) pdffile->Get("sk_RHC1Re_numubarcccoh");
  if(!sk_1Re_numubarccoth) sk_1Re_numubarccoth = (TH1D*) pdffile->Get("sk_RHC1Re_numubarccoth");
  if(!sk_1Re_numubarnc1pi0) sk_1Re_numubarnc1pi0 = (TH1D*) pdffile->Get("sk_RHC1Re_numubarnc1pi0");
  if(!sk_1Re_numubarnccoh) sk_1Re_numubarnccoh = (TH1D*) pdffile->Get("sk_RHC1Re_numubarnccoh");
  if(!sk_1Re_numubarncoth) sk_1Re_numubarncoth = (TH1D*) pdffile->Get("sk_RHC1Re_numubarncoth");
  if(!sk_1Re_nueccqe) sk_1Re_nueccqe = (TH1D*) pdffile->Get("sk_RHC1Re_nueccqe");
  if(!sk_1Re_nuecc1pi) sk_1Re_nuecc1pi = (TH1D*) pdffile->Get("sk_RHC1Re_nuecc1pi");
  if(!sk_1Re_nuecccoh) sk_1Re_nuecccoh = (TH1D*) pdffile->Get("sk_RHC1Re_nuecccoh");
  if(!sk_1Re_nueccoth) sk_1Re_nueccoth = (TH1D*) pdffile->Get("sk_RHC1Re_nueccoth");
  if(!sk_1Re_nuenc1pi0) sk_1Re_nuenc1pi0 = (TH1D*) pdffile->Get("sk_RHC1Re_nuenc1pi0");
  if(!sk_1Re_nuenccoh) sk_1Re_nuenccoh = (TH1D*) pdffile->Get("sk_RHC1Re_nuenccoh");
  if(!sk_1Re_nuencoth) sk_1Re_nuencoth = (TH1D*) pdffile->Get("sk_RHC1Re_nuencoth");
  if(!sk_1Re_nuebarccqe) sk_1Re_nuebarccqe = (TH1D*) pdffile->Get("sk_RHC1Re_nuebarccqe");
  if(!sk_1Re_nuebarcc1pi) sk_1Re_nuebarcc1pi = (TH1D*) pdffile->Get("sk_RHC1Re_nuebarcc1pi");
  if(!sk_1Re_nuebarcccoh) sk_1Re_nuebarcccoh = (TH1D*) pdffile->Get("sk_RHC1Re_nuebarcccoh");
  if(!sk_1Re_nuebarccoth) sk_1Re_nuebarccoth = (TH1D*) pdffile->Get("sk_RHC1Re_nuebarccoth");
  if(!sk_1Re_nuebarnc1pi0) sk_1Re_nuebarnc1pi0 = (TH1D*) pdffile->Get("sk_RHC1Re_nuebarnc1pi0");
  if(!sk_1Re_nuebarnccoh) sk_1Re_nuebarnccoh = (TH1D*) pdffile->Get("sk_RHC1Re_nuebarnccoh");
  if(!sk_1Re_nuebarncoth) sk_1Re_nuebarncoth = (TH1D*) pdffile->Get("sk_RHC1Re_nuebarncoth");
  if(!sk_1Re_oscnumubarccqe) sk_1Re_oscnumubarccqe = (TH1D*) pdffile->Get("sk_RHC1Re_oscnumubarccqe");
  if(!sk_1Re_oscnumubarcc1pi) sk_1Re_oscnumubarcc1pi = (TH1D*) pdffile->Get("sk_RHC1Re_oscnumubarcc1pi");
  if(!sk_1Re_oscnumubarcccoh) sk_1Re_oscnumubarcccoh = (TH1D*) pdffile->Get("sk_RHC1Re_oscnumubarcccoh");
  if(!sk_1Re_oscnumubarccoth) sk_1Re_oscnumubarccoth = (TH1D*) pdffile->Get("sk_RHC1Re_oscnumubarccoth");
  if(!sk_1Re_oscnuebarccqe) sk_1Re_oscnuebarccqe = (TH1D*) pdffile->Get("sk_RHC1Re_oscnuebarccqe");
  if(!sk_1Re_oscnuebarcc1pi) sk_1Re_oscnuebarcc1pi = (TH1D*) pdffile->Get("sk_RHC1Re_oscnuebarcc1pi");
  if(!sk_1Re_oscnuebarcccoh) sk_1Re_oscnuebarcccoh = (TH1D*) pdffile->Get("sk_RHC1Re_oscnuebarcccoh");
  if(!sk_1Re_oscnuebarccoth) sk_1Re_oscnuebarccoth = (TH1D*) pdffile->Get("sk_RHC1Re_oscnuebarccoth");
  if(!sk_1Re_oscnueccqe) sk_1Re_oscnueccqe = (TH1D*) pdffile->Get("sk_RHC1Re_oscnueccqe");
  if(!sk_1Re_oscnuecc1pi) sk_1Re_oscnuecc1pi = (TH1D*) pdffile->Get("sk_RHC1Re_oscnuecc1pi");
  if(!sk_1Re_oscnuecccoh) sk_1Re_oscnuecccoh = (TH1D*) pdffile->Get("sk_RHC1Re_oscnuecccoh");
  if(!sk_1Re_oscnueccoth) sk_1Re_oscnueccoth = (TH1D*) pdffile->Get("sk_RHC1Re_oscnueccoth");

  //
  // Clone pdfs to have local versions. Rebin if requested.
  //

  TH1D * h_sk_1Re_total         = 0;
  TH1D * h_sk_1Re_numuccqe      = 0;
  TH1D * h_sk_1Re_numucc1pi     = 0;
  TH1D * h_sk_1Re_numucccoh     = 0;
  TH1D * h_sk_1Re_numuccoth     = 0;
  TH1D * h_sk_1Re_numunc1pi0    = 0;
  TH1D * h_sk_1Re_numunccoh    = 0;
  TH1D * h_sk_1Re_numuncoth     = 0;
  TH1D * h_sk_1Re_oscnumuccqe   = 0;
  TH1D * h_sk_1Re_oscnumucc1pi  = 0;
  TH1D * h_sk_1Re_oscnumucccoh  = 0;
  TH1D * h_sk_1Re_oscnumuccoth  = 0;
  TH1D * h_sk_1Re_numubarccqe   = 0;
  TH1D * h_sk_1Re_numubarcc1pi  = 0;
  TH1D * h_sk_1Re_numubarcccoh  = 0;
  TH1D * h_sk_1Re_numubarccoth  = 0;
  TH1D * h_sk_1Re_numubarnc1pi0 = 0;
  TH1D * h_sk_1Re_numubarnccoh = 0;
  TH1D * h_sk_1Re_numubarncoth  = 0;
  TH1D * h_sk_1Re_nueccqe	 = 0;
  TH1D * h_sk_1Re_nuecc1pi      = 0;
  TH1D * h_sk_1Re_nuecccoh      = 0;
  TH1D * h_sk_1Re_nueccoth      = 0;
  TH1D * h_sk_1Re_nuenc1pi0     = 0;
  TH1D * h_sk_1Re_nuenccoh     = 0;
  TH1D * h_sk_1Re_nuencoth      = 0;
  TH1D * h_sk_1Re_nuebarccqe    = 0;
  TH1D * h_sk_1Re_nuebarcc1pi   = 0;
  TH1D * h_sk_1Re_nuebarcccoh   = 0;
  TH1D * h_sk_1Re_nuebarccoth   = 0;
  TH1D * h_sk_1Re_nuebarnc1pi0  = 0;
  TH1D * h_sk_1Re_nuebarnccoh  = 0;
  TH1D * h_sk_1Re_nuebarncoth   = 0;
  TH1D * h_sk_1Re_oscnumubarccqe    = 0;
  TH1D * h_sk_1Re_oscnumubarcc1pi   = 0;
  TH1D * h_sk_1Re_oscnumubarcccoh   = 0;
  TH1D * h_sk_1Re_oscnumubarccoth   = 0;
  TH1D * h_sk_1Re_oscnueccqe    = 0;
  TH1D * h_sk_1Re_oscnuecc1pi   = 0;
  TH1D * h_sk_1Re_oscnuecccoh   = 0;
  TH1D * h_sk_1Re_oscnueccoth   = 0;
  TH1D * h_sk_1Re_oscnuebarccqe    = 0;
  TH1D * h_sk_1Re_oscnuebarcc1pi   = 0;
  TH1D * h_sk_1Re_oscnuebarcccoh   = 0;
  TH1D * h_sk_1Re_oscnuebarccoth   = 0;

  h_sk_1Re_total         = new TH1D(*sk_1Retotal);
  h_sk_1Re_numuccqe      = new TH1D(*sk_1Re_numuccqe);
  h_sk_1Re_numucc1pi     = new TH1D(*sk_1Re_numucc1pi);
  h_sk_1Re_numucccoh     = new TH1D(*sk_1Re_numucccoh);
  h_sk_1Re_numuccoth     = new TH1D(*sk_1Re_numuccoth);
  h_sk_1Re_numunc1pi0    = new TH1D(*sk_1Re_numunc1pi0);
  h_sk_1Re_numunccoh    = new TH1D(*sk_1Re_numunccoh);
  h_sk_1Re_numuncoth     = new TH1D(*sk_1Re_numuncoth);
  h_sk_1Re_oscnumuccqe      = new TH1D(*sk_1Re_oscnumuccqe);
  h_sk_1Re_oscnumucc1pi     = new TH1D(*sk_1Re_oscnumucc1pi);
  h_sk_1Re_oscnumucccoh     = new TH1D(*sk_1Re_oscnumucccoh);
  h_sk_1Re_oscnumuccoth     = new TH1D(*sk_1Re_oscnumuccoth);
  h_sk_1Re_numubarccqe   = new TH1D(*sk_1Re_numubarccqe);
  h_sk_1Re_numubarcc1pi  = new TH1D(*sk_1Re_numubarcc1pi);
  h_sk_1Re_numubarcccoh  = new TH1D(*sk_1Re_numubarcccoh);
  h_sk_1Re_numubarccoth  = new TH1D(*sk_1Re_numubarccoth);
  h_sk_1Re_numubarnc1pi0 = new TH1D(*sk_1Re_numubarnc1pi0);
  h_sk_1Re_numubarnccoh = new TH1D(*sk_1Re_numubarnccoh);
  h_sk_1Re_numubarncoth  = new TH1D(*sk_1Re_numubarncoth);
  h_sk_1Re_nueccqe	     = new TH1D(*sk_1Re_nueccqe);	
  h_sk_1Re_nuecc1pi      = new TH1D(*sk_1Re_nuecc1pi);
  h_sk_1Re_nuecccoh      = new TH1D(*sk_1Re_nuecccoh);
  h_sk_1Re_nueccoth      = new TH1D(*sk_1Re_nueccoth);
  h_sk_1Re_nuenc1pi0     = new TH1D(*sk_1Re_nuenc1pi0);
  h_sk_1Re_nuenccoh     = new TH1D(*sk_1Re_nuenccoh);
  h_sk_1Re_nuencoth      = new TH1D(*sk_1Re_nuencoth);
  h_sk_1Re_nuebarccqe    = new TH1D(*sk_1Re_nuebarccqe);
  h_sk_1Re_nuebarcc1pi   = new TH1D(*sk_1Re_nuebarcc1pi);
  h_sk_1Re_nuebarcccoh   = new TH1D(*sk_1Re_nuebarcccoh);
  h_sk_1Re_nuebarccoth   = new TH1D(*sk_1Re_nuebarccoth);
  h_sk_1Re_nuebarnc1pi0  = new TH1D(*sk_1Re_nuebarnc1pi0);
  h_sk_1Re_nuebarnccoh  = new TH1D(*sk_1Re_nuebarnccoh);
  h_sk_1Re_nuebarncoth   = new TH1D(*sk_1Re_nuebarncoth);
  h_sk_1Re_oscnumubarccqe    = new TH1D(*sk_1Re_oscnumubarccqe);
  h_sk_1Re_oscnumubarcc1pi   = new TH1D(*sk_1Re_oscnumubarcc1pi);
  h_sk_1Re_oscnumubarcccoh   = new TH1D(*sk_1Re_oscnumubarcccoh);
  h_sk_1Re_oscnumubarccoth   = new TH1D(*sk_1Re_oscnumubarccoth);
  h_sk_1Re_oscnueccqe    = new TH1D(*sk_1Re_oscnueccqe);
  h_sk_1Re_oscnuecc1pi   = new TH1D(*sk_1Re_oscnuecc1pi);
  h_sk_1Re_oscnuecccoh   = new TH1D(*sk_1Re_oscnuecccoh);
  h_sk_1Re_oscnueccoth   = new TH1D(*sk_1Re_oscnueccoth);
  h_sk_1Re_oscnuebarccqe    = new TH1D(*sk_1Re_oscnuebarccqe);
  h_sk_1Re_oscnuebarcc1pi   = new TH1D(*sk_1Re_oscnuebarcc1pi);
  h_sk_1Re_oscnuebarcccoh   = new TH1D(*sk_1Re_oscnuebarcccoh);
  h_sk_1Re_oscnuebarccoth   = new TH1D(*sk_1Re_oscnuebarccoth);

//  pdffile->Close();

  //
  // Now re-group pdfs in fewer categories
  //
  TH1D * h_mc_all           = 0;
  TH1D * h_mc_oscnueccqe    = 0;
  TH1D * h_mc_oscnueccnonqe = 0;
  TH1D * h_mc_nuecc         = 0;
  TH1D * h_mc_numucc        = 0;
  TH1D * h_mc_nc            = 0;

  h_mc_all = new TH1D(*h_sk_1Re_total);

  h_mc_oscnueccqe  = new TH1D(*h_sk_1Re_oscnueccqe);
  h_mc_oscnueccqe -> Add (h_sk_1Re_oscnuebarccqe);

  h_mc_oscnueccnonqe = new TH1D(*h_sk_1Re_oscnuecc1pi);
  h_mc_oscnueccnonqe -> Add (h_sk_1Re_oscnuecccoh);
  h_mc_oscnueccnonqe -> Add (h_sk_1Re_oscnueccoth);
  h_mc_oscnueccnonqe -> Add (h_sk_1Re_oscnuebarcc1pi);
  h_mc_oscnueccnonqe -> Add (h_sk_1Re_oscnuebarcccoh);
  h_mc_oscnueccnonqe -> Add (h_sk_1Re_oscnuebarccoth);

  h_mc_nuecc = new TH1D(*h_sk_1Re_nueccqe);
  h_mc_nuecc -> Add (h_sk_1Re_nuecc1pi);
  h_mc_nuecc -> Add (h_sk_1Re_nuecccoh);
  h_mc_nuecc -> Add (h_sk_1Re_nueccoth);
  h_mc_nuecc -> Add (h_sk_1Re_nuebarccqe);
  h_mc_nuecc -> Add (h_sk_1Re_nuebarcc1pi);
  h_mc_nuecc -> Add (h_sk_1Re_nuebarcccoh);
  h_mc_nuecc -> Add (h_sk_1Re_nuebarccoth);
  
  h_mc_numucc = new TH1D(*h_sk_1Re_numuccqe);
  h_mc_numucc -> Add (h_sk_1Re_numucc1pi);
  h_mc_numucc -> Add (h_sk_1Re_numucccoh);
  h_mc_numucc -> Add (h_sk_1Re_numuccoth);
  h_mc_numucc -> Add (h_sk_1Re_oscnumuccqe);
  h_mc_numucc -> Add (h_sk_1Re_oscnumucc1pi);
  h_mc_numucc -> Add (h_sk_1Re_oscnumucccoh);
  h_mc_numucc -> Add (h_sk_1Re_oscnumuccoth);
  h_mc_numucc -> Add (h_sk_1Re_numubarccqe);
  h_mc_numucc -> Add (h_sk_1Re_numubarcc1pi);
  h_mc_numucc -> Add (h_sk_1Re_numubarcccoh);
  h_mc_numucc -> Add (h_sk_1Re_numubarccoth);
  h_mc_numucc -> Add (h_sk_1Re_oscnumubarccqe);
  h_mc_numucc -> Add (h_sk_1Re_oscnumubarcc1pi);
  h_mc_numucc -> Add (h_sk_1Re_oscnumubarcccoh);
  h_mc_numucc -> Add (h_sk_1Re_oscnumubarccoth);

  h_mc_nc = new TH1D(*h_sk_1Re_numunc1pi0);
  h_mc_nc     -> Add (h_sk_1Re_numunccoh);
  h_mc_nc     -> Add (h_sk_1Re_numuncoth);
  h_mc_nc     -> Add (h_sk_1Re_numubarnc1pi0);
  h_mc_nc     -> Add (h_sk_1Re_numubarnccoh);
  h_mc_nc     -> Add (h_sk_1Re_numubarncoth);
  h_mc_nc     -> Add (h_sk_1Re_nuenc1pi0);
  h_mc_nc     -> Add (h_sk_1Re_nuenccoh);
  h_mc_nc     -> Add (h_sk_1Re_nuencoth);
  h_mc_nc     -> Add (h_sk_1Re_nuebarnc1pi0);
  h_mc_nc     -> Add (h_sk_1Re_nuebarnccoh);
  h_mc_nc     -> Add (h_sk_1Re_nuebarncoth);

  //
  // If a second pdffile was specified, get the total only
  //

  TH1D * h_mc_all_alt = 0;
  if(altpdffilename) {
    // get pdf
    TFile * altpdffile = new TFile(altpdffilename, "read");
    altpdffile->cd();
    TH1D * sk_1Retotal_alt = (TH1D*) pdffile->Get("sk_1Retotal");
    // create local copy, rebinned or not
    if(rebin) {
       h_mc_all_alt = (TH1D*) sk_1Retotal_alt -> Rebin(kNErecoBins1Re, "h_sk_1Re_total_alt", kErecoBinEdges1Re);
    } else {
       h_mc_all_alt = new TH1D(*sk_1Retotal_alt);
    }
    altpdffile->Close();
  }

  //
  // If a datafile was specified, get spectrum
  //

  TH1D * h_data = 0;

  if(datafilename) {
     TFile * datafile = new TFile(datafilename, "read");
     datafile->cd();
     TH1D * hd = (TH1D *) datafile->Get("obs_spectrum");
     if(rebin) {
         h_data = (TH1D*) hd -> Rebin(kNErecoBins1Re, "histo_data", kErecoBinEdges1Re);
     } else {
         h_data = new TH1D (*hd);
     }
     cout << "Read " << h_data->GetEntries() << " events " << endl;
  }


  //
  // Format
  //

  h_mc_all         -> SetLineWidth(2);
  h_mc_all         -> SetLineColor(kBlack);
  h_mc_all         -> SetLineStyle(kSolid);
  h_mc_oscnueccqe -> SetLineColor(kBlack);
  h_mc_oscnueccqe -> SetFillColor(kWhite);
  h_mc_oscnueccqe -> SetFillStyle(3354);
  h_mc_oscnueccqe -> SetLineColor(kGray+3);
  h_mc_oscnueccqe -> SetFillColor(kGray+3);
  h_mc_oscnueccnonqe   -> SetLineColor(kGreen+1);
  h_mc_oscnueccnonqe   -> SetFillColor(kGreen+1);
  h_mc_nuecc   -> SetLineColor(kCyan+1);
  h_mc_nuecc   -> SetFillColor(kCyan+1);
  h_mc_numucc -> SetLineColor(kAzure-3);
  h_mc_numucc -> SetFillColor(kAzure-3);
  h_mc_nc          -> SetLineColor(kRed);
  h_mc_nc          -> SetFillColor(kRed);
  if(h_mc_all_alt) {
     h_mc_all_alt -> SetLineWidth(2);
     h_mc_all_alt -> SetLineColor(kRed);
     h_mc_all_alt -> SetLineStyle(kSolid);
  }
  if(h_data) {
     h_data -> SetMarkerStyle(8);
     h_data -> SetMarkerSize(1);
     h_data -> SetLineStyle(kDashed);
     h_data -> SetLineWidth(1);
  }

  // turn into events per 50 MeV
  for(int ix = 1; ix <= h_mc_all->GetNbinsX(); ix++) {
    h_mc_all->SetBinContent(ix, (h_mc_all->GetBinContent(ix) / h_mc_all->GetBinWidth(ix)) * 0.05);
    h_mc_oscnueccqe->SetBinContent(ix, (h_mc_oscnueccqe->GetBinContent(ix) / h_mc_oscnueccqe->GetBinWidth(ix)) * 0.05);
    h_mc_oscnueccnonqe->SetBinContent(ix, (h_mc_oscnueccnonqe->GetBinContent(ix) / h_mc_oscnueccnonqe->GetBinWidth(ix)) * 0.05);
    h_mc_nuecc->SetBinContent(ix, (h_mc_nuecc->GetBinContent(ix) / h_mc_nuecc->GetBinWidth(ix)) * 0.05);
    h_mc_numucc->SetBinContent(ix, (h_mc_numucc->GetBinContent(ix) / h_mc_numucc->GetBinWidth(ix)) * 0.05);
    if(h_mc_all_alt)
      h_mc_all_alt->SetBinContent(ix, (h_mc_all_alt->GetBinContent(ix) / h_mc_all_alt->GetBinWidth(ix)) * 0.05);
    if(h_data)
      h_data->SetBinContent(ix, (h_data->GetBinContent(ix) / h_data->GetBinWidth(ix)) * 0.05);
  }//ix

  // stack together components of the first pdf
  THStack * stack = new THStack("stack","spectrum contributions");
  stack->Add(h_mc_nc);
  stack->Add(h_mc_numucc);
  stack->Add(h_mc_nuecc);
  stack->Add(h_mc_oscnueccnonqe);
  stack->Add(h_mc_oscnueccqe);

  //
  // Start plotting
  //

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  c->SetLeftMargin(0.15);
  TH1F * hframe = (TH1F*)c->DrawFrame(0,0,emax,nmax);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitleOffset(1.2);
  hframe->GetYaxis()->SetTitle("Events per 0.050 GeV");

  h_mc_all->Draw("same");
  if(show_components) {
     stack->Draw("same");
  }

  if(h_mc_all_alt) {
     h_mc_all_alt->Draw("same");
  }

  if(h_data) {
     h_data->Draw("PE03X0 same");
  }

  // 
  // Add legend
  // 
  TLegend * legend = new TLegend(0.5,0.6,0.85,0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->SetTextFont(132);
  if(print_header)
    legend->SetHeader("1 #mu-like ring");
  if(h_data) {
     legend->AddEntry(h_data,     Form("%s",datalabel),            "LP");
  }
  //legend->AddEntry(h_mc_all, Form("%s (total)",pdflabel),"L");
  if(show_components) {
     legend->AddEntry(h_mc_oscnueccqe,  Form("%s Osc. #nu_{e} CCQE",pdflabel),         "LF");
     legend->AddEntry(h_mc_oscnueccnonqe,    Form("%s Osc. #nu_{e} CCnonQE",pdflabel),        "LF");
     legend->AddEntry(h_mc_nuecc,    Form("%s Bkg. #nu_{e}+#bar{#nu_{e}}CC",pdflabel),        "LF");
     legend->AddEntry(h_mc_numucc,  Form("%s #nu_{#mu}+#bar{#nu_{#mu}} CC",pdflabel), "LF");
     legend->AddEntry(h_mc_nc,           Form("%s NC",pdflabel),                        "LF");
  }
  legend->Draw();

  c->SaveAs("spectrum.eps");
  c->SaveAs("spectrum.pdf");
}

