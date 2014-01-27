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

void plot_1ringmu_spectrum_rebin_regroup
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
  double emax = 6.0,      

  /* max of N_{events} axis in plot  */
  double nmax = 20.0,

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
  TH1D * sk_1Rmutotal          = (TH1D*) pdffile->Get("sk_1Rmutotal");
  TH1D * sk_1Rmu_numuccqe      = (TH1D*) pdffile->Get("sk_1Rmu_numuccqe");
  TH1D * sk_1Rmu_numucc1pi     = (TH1D*) pdffile->Get("sk_1Rmu_numucc1pi");
  TH1D * sk_1Rmu_numucccoh     = (TH1D*) pdffile->Get("sk_1Rmu_numucccoh");
  TH1D * sk_1Rmu_numuccoth     = (TH1D*) pdffile->Get("sk_1Rmu_numuccoth");
  TH1D * sk_1Rmu_numunc1piC    = (TH1D*) pdffile->Get("sk_1Rmu_numunc1piC");
  TH1D * sk_1Rmu_numuncoth     = (TH1D*) pdffile->Get("sk_1Rmu_numuncoth");
  TH1D * sk_1Rmu_oscnumuccqe   = (TH1D*) pdffile->Get("sk_1Rmu_oscnumuccqe");
  TH1D * sk_1Rmu_oscnumucc1pi  = (TH1D*) pdffile->Get("sk_1Rmu_oscnumucc1pi");
  TH1D * sk_1Rmu_oscnumucccoh  = (TH1D*) pdffile->Get("sk_1Rmu_oscnumucccoh");
  TH1D * sk_1Rmu_oscnumuccoth  = (TH1D*) pdffile->Get("sk_1Rmu_oscnumuccoth");
  TH1D * sk_1Rmu_numubarccqe   = (TH1D*) pdffile->Get("sk_1Rmu_numubarccqe");
  TH1D * sk_1Rmu_numubarcc1pi  = (TH1D*) pdffile->Get("sk_1Rmu_numubarcc1pi");
  TH1D * sk_1Rmu_numubarcccoh  = (TH1D*) pdffile->Get("sk_1Rmu_numubarcccoh");
  TH1D * sk_1Rmu_numubarccoth  = (TH1D*) pdffile->Get("sk_1Rmu_numubarccoth");
  TH1D * sk_1Rmu_numubarnc1piC = (TH1D*) pdffile->Get("sk_1Rmu_numubarnc1piC");
  TH1D * sk_1Rmu_numubarncoth  = (TH1D*) pdffile->Get("sk_1Rmu_numubarncoth");
  TH1D * sk_1Rmu_nueccqe       = (TH1D*) pdffile->Get("sk_1Rmu_nueccqe");
  TH1D * sk_1Rmu_nuecc1pi      = (TH1D*) pdffile->Get("sk_1Rmu_nuecc1pi");
  TH1D * sk_1Rmu_nuecccoh      = (TH1D*) pdffile->Get("sk_1Rmu_nuecccoh");
  TH1D * sk_1Rmu_nueccoth      = (TH1D*) pdffile->Get("sk_1Rmu_nueccoth");
  TH1D * sk_1Rmu_nuenc1piC     = (TH1D*) pdffile->Get("sk_1Rmu_nuenc1piC");
  TH1D * sk_1Rmu_nuencoth      = (TH1D*) pdffile->Get("sk_1Rmu_nuencoth");
  TH1D * sk_1Rmu_nuebarccqe    = (TH1D*) pdffile->Get("sk_1Rmu_nuebarccqe");
  TH1D * sk_1Rmu_nuebarcc1pi   = (TH1D*) pdffile->Get("sk_1Rmu_nuebarcc1pi");
  TH1D * sk_1Rmu_nuebarcccoh   = (TH1D*) pdffile->Get("sk_1Rmu_nuebarcccoh");
  TH1D * sk_1Rmu_nuebarccoth   = (TH1D*) pdffile->Get("sk_1Rmu_nuebarccoth");
  TH1D * sk_1Rmu_nuebarnc1piC  = (TH1D*) pdffile->Get("sk_1Rmu_nuebarnc1piC");
  TH1D * sk_1Rmu_nuebarncoth   = (TH1D*) pdffile->Get("sk_1Rmu_nuebarncoth");
  TH1D * sk_1Rmu_oscnumubarccqe    = (TH1D*) pdffile->Get("sk_1Rmu_oscnumubarccqe");
  TH1D * sk_1Rmu_oscnumubarcc1pi   = (TH1D*) pdffile->Get("sk_1Rmu_oscnumubarcc1pi");
  TH1D * sk_1Rmu_oscnumubarcccoh   = (TH1D*) pdffile->Get("sk_1Rmu_oscnumubarcccoh");
  TH1D * sk_1Rmu_oscnumubarccoth   = (TH1D*) pdffile->Get("sk_1Rmu_oscnumubarccoth");
  TH1D * sk_1Rmu_oscnueccqe    = (TH1D*) pdffile->Get("sk_1Rmu_oscnueccqe");
  TH1D * sk_1Rmu_oscnuecc1pi   = (TH1D*) pdffile->Get("sk_1Rmu_oscnuecc1pi");
  TH1D * sk_1Rmu_oscnuecccoh   = (TH1D*) pdffile->Get("sk_1Rmu_oscnuecccoh");
  TH1D * sk_1Rmu_oscnueccoth   = (TH1D*) pdffile->Get("sk_1Rmu_oscnueccoth");
  TH1D * sk_1Rmu_oscnuebarccqe    = (TH1D*) pdffile->Get("sk_1Rmu_oscnuebarccqe");
  TH1D * sk_1Rmu_oscnuebarcc1pi   = (TH1D*) pdffile->Get("sk_1Rmu_oscnuebarcc1pi");
  TH1D * sk_1Rmu_oscnuebarcccoh   = (TH1D*) pdffile->Get("sk_1Rmu_oscnuebarcccoh");
  TH1D * sk_1Rmu_oscnuebarccoth   = (TH1D*) pdffile->Get("sk_1Rmu_oscnuebarccoth");

  if(!sk_1Rmutotal) sk_1Rmutotal = (TH1D*) pdffile->Get("sk_RHC1Rmutotal");
  if(!sk_1Rmu_numuccqe) sk_1Rmu_numuccqe = (TH1D*) pdffile->Get("sk_RHC1Rmu_numuccqe");
  if(!sk_1Rmu_numucc1pi) sk_1Rmu_numucc1pi = (TH1D*) pdffile->Get("sk_RHC1Rmu_numucc1pi");
  if(!sk_1Rmu_numucccoh) sk_1Rmu_numucccoh = (TH1D*) pdffile->Get("sk_RHC1Rmu_numucccoh");
  if(!sk_1Rmu_numuccoth) sk_1Rmu_numuccoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_numuccoth");
  if(!sk_1Rmu_numunc1piC) sk_1Rmu_numunc1piC = (TH1D*) pdffile->Get("sk_RHC1Rmu_numunc1piC");
  if(!sk_1Rmu_numuncoth) sk_1Rmu_numuncoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_numuncoth");
  if(!sk_1Rmu_oscnumuccqe) sk_1Rmu_oscnumuccqe = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnumuccqe");
  if(!sk_1Rmu_oscnumucc1pi) sk_1Rmu_oscnumucc1pi = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnumucc1pi");
  if(!sk_1Rmu_oscnumucccoh) sk_1Rmu_oscnumucccoh = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnumucccoh");
  if(!sk_1Rmu_oscnumuccoth) sk_1Rmu_oscnumuccoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnumuccoth");
  if(!sk_1Rmu_numubarccqe) sk_1Rmu_numubarccqe = (TH1D*) pdffile->Get("sk_RHC1Rmu_numubarccqe");
  if(!sk_1Rmu_numubarcc1pi) sk_1Rmu_numubarcc1pi = (TH1D*) pdffile->Get("sk_RHC1Rmu_numubarcc1pi");
  if(!sk_1Rmu_numubarcccoh) sk_1Rmu_numubarcccoh = (TH1D*) pdffile->Get("sk_RHC1Rmu_numubarcccoh");
  if(!sk_1Rmu_numubarccoth) sk_1Rmu_numubarccoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_numubarccoth");
  if(!sk_1Rmu_numubarnc1piC) sk_1Rmu_numubarnc1piC = (TH1D*) pdffile->Get("sk_RHC1Rmu_numubarnc1piC");
  if(!sk_1Rmu_numubarncoth) sk_1Rmu_numubarncoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_numubarncoth");
  if(!sk_1Rmu_nueccqe) sk_1Rmu_nueccqe = (TH1D*) pdffile->Get("sk_RHC1Rmu_nueccqe");
  if(!sk_1Rmu_nuecc1pi) sk_1Rmu_nuecc1pi = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuecc1pi");
  if(!sk_1Rmu_nuecccoh) sk_1Rmu_nuecccoh = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuecccoh");
  if(!sk_1Rmu_nueccoth) sk_1Rmu_nueccoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_nueccoth");
  if(!sk_1Rmu_nuenc1piC) sk_1Rmu_nuenc1piC = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuenc1piC");
  if(!sk_1Rmu_nuencoth) sk_1Rmu_nuencoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuencoth");
  if(!sk_1Rmu_nuebarccqe) sk_1Rmu_nuebarccqe = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuebarccqe");
  if(!sk_1Rmu_nuebarcc1pi) sk_1Rmu_nuebarcc1pi = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuebarcc1pi");
  if(!sk_1Rmu_nuebarcccoh) sk_1Rmu_nuebarcccoh = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuebarcccoh");
  if(!sk_1Rmu_nuebarccoth) sk_1Rmu_nuebarccoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuebarccoth");
  if(!sk_1Rmu_nuebarnc1piC) sk_1Rmu_nuebarnc1piC = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuebarnc1piC");
  if(!sk_1Rmu_nuebarncoth) sk_1Rmu_nuebarncoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_nuebarncoth");
  if(!sk_1Rmu_oscnumubarccqe) sk_1Rmu_oscnumubarccqe = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnumubarccqe");
  if(!sk_1Rmu_oscnumubarcc1pi) sk_1Rmu_oscnumubarcc1pi = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnumubarcc1pi");
  if(!sk_1Rmu_oscnumubarcccoh) sk_1Rmu_oscnumubarcccoh = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnumubarcccoh");
  if(!sk_1Rmu_oscnumubarccoth) sk_1Rmu_oscnumubarccoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnumubarccoth");
  if(!sk_1Rmu_oscnuebarccqe) sk_1Rmu_oscnuebarccqe = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnuebarccqe");
  if(!sk_1Rmu_oscnuebarcc1pi) sk_1Rmu_oscnuebarcc1pi = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnuebarcc1pi");
  if(!sk_1Rmu_oscnuebarcccoh) sk_1Rmu_oscnuebarcccoh = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnuebarcccoh");
  if(!sk_1Rmu_oscnuebarccoth) sk_1Rmu_oscnuebarccoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnuebarccoth");
  if(!sk_1Rmu_oscnueccqe) sk_1Rmu_oscnueccqe = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnueccqe");
  if(!sk_1Rmu_oscnuecc1pi) sk_1Rmu_oscnuecc1pi = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnuecc1pi");
  if(!sk_1Rmu_oscnuecccoh) sk_1Rmu_oscnuecccoh = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnuecccoh");
  if(!sk_1Rmu_oscnueccoth) sk_1Rmu_oscnueccoth = (TH1D*) pdffile->Get("sk_RHC1Rmu_oscnueccoth");

  //
  // Clone pdfs to have local versions. Rebin if requested.
  //

  TH1D * h_sk_1Rmu_total         = 0;
  TH1D * h_sk_1Rmu_numuccqe      = 0;
  TH1D * h_sk_1Rmu_numucc1pi     = 0;
  TH1D * h_sk_1Rmu_numucccoh     = 0;
  TH1D * h_sk_1Rmu_numuccoth     = 0;
  TH1D * h_sk_1Rmu_numunc1piC    = 0;
  TH1D * h_sk_1Rmu_numuncoth     = 0;
  TH1D * h_sk_1Rmu_oscnumuccqe   = 0;
  TH1D * h_sk_1Rmu_oscnumucc1pi  = 0;
  TH1D * h_sk_1Rmu_oscnumucccoh  = 0;
  TH1D * h_sk_1Rmu_oscnumuccoth  = 0;
  TH1D * h_sk_1Rmu_numubarccqe   = 0;
  TH1D * h_sk_1Rmu_numubarcc1pi  = 0;
  TH1D * h_sk_1Rmu_numubarcccoh  = 0;
  TH1D * h_sk_1Rmu_numubarccoth  = 0;
  TH1D * h_sk_1Rmu_numubarnc1piC = 0;
  TH1D * h_sk_1Rmu_numubarncoth  = 0;
  TH1D * h_sk_1Rmu_nueccqe	 = 0;
  TH1D * h_sk_1Rmu_nuecc1pi      = 0;
  TH1D * h_sk_1Rmu_nuecccoh      = 0;
  TH1D * h_sk_1Rmu_nueccoth      = 0;
  TH1D * h_sk_1Rmu_nuenc1piC     = 0;
  TH1D * h_sk_1Rmu_nuencoth      = 0;
  TH1D * h_sk_1Rmu_nuebarccqe    = 0;
  TH1D * h_sk_1Rmu_nuebarcc1pi   = 0;
  TH1D * h_sk_1Rmu_nuebarcccoh   = 0;
  TH1D * h_sk_1Rmu_nuebarccoth   = 0;
  TH1D * h_sk_1Rmu_nuebarnc1piC  = 0;
  TH1D * h_sk_1Rmu_nuebarncoth   = 0;
  TH1D * h_sk_1Rmu_oscnumubarccqe    = 0;
  TH1D * h_sk_1Rmu_oscnumubarcc1pi   = 0;
  TH1D * h_sk_1Rmu_oscnumubarcccoh   = 0;
  TH1D * h_sk_1Rmu_oscnumubarccoth   = 0;
  TH1D * h_sk_1Rmu_oscnueccqe    = 0;
  TH1D * h_sk_1Rmu_oscnuecc1pi   = 0;
  TH1D * h_sk_1Rmu_oscnuecccoh   = 0;
  TH1D * h_sk_1Rmu_oscnueccoth   = 0;
  TH1D * h_sk_1Rmu_oscnuebarccqe    = 0;
  TH1D * h_sk_1Rmu_oscnuebarcc1pi   = 0;
  TH1D * h_sk_1Rmu_oscnuebarcccoh   = 0;
  TH1D * h_sk_1Rmu_oscnuebarccoth   = 0;

  if(rebin) {
     h_sk_1Rmu_total         = (TH1D*) sk_1Rmutotal          -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_total",         kErecoBinEdges1Rmu);
     h_sk_1Rmu_numuccqe      = (TH1D*) sk_1Rmu_numuccqe      -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numuccqe",      kErecoBinEdges1Rmu);
     h_sk_1Rmu_numucc1pi     = (TH1D*) sk_1Rmu_numucc1pi     -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numucc1pi",     kErecoBinEdges1Rmu);
     h_sk_1Rmu_numucccoh     = (TH1D*) sk_1Rmu_numucccoh     -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numucccoh",     kErecoBinEdges1Rmu);
     h_sk_1Rmu_numuccoth     = (TH1D*) sk_1Rmu_numuccoth     -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numuccoth",     kErecoBinEdges1Rmu);
     h_sk_1Rmu_numunc1piC    = (TH1D*) sk_1Rmu_numunc1piC    -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numunc1piC",    kErecoBinEdges1Rmu);
     h_sk_1Rmu_numuncoth     = (TH1D*) sk_1Rmu_numuncoth     -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numuncoth",     kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnumuccqe      = (TH1D*) sk_1Rmu_oscnumuccqe      -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnumuccqe",      kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnumucc1pi     = (TH1D*) sk_1Rmu_oscnumucc1pi     -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnumucc1pi",     kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnumucccoh     = (TH1D*) sk_1Rmu_oscnumucccoh     -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnumucccoh",     kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnumuccoth     = (TH1D*) sk_1Rmu_oscnumuccoth     -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnumuccoth",     kErecoBinEdges1Rmu);
     h_sk_1Rmu_numubarccqe   = (TH1D*) sk_1Rmu_numubarccqe   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numubarccqe",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_numubarcc1pi  = (TH1D*) sk_1Rmu_numubarcc1pi  -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numubarcc1pi",  kErecoBinEdges1Rmu);
     h_sk_1Rmu_numubarcccoh  = (TH1D*) sk_1Rmu_numubarcccoh  -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numubarcccoh",  kErecoBinEdges1Rmu);
     h_sk_1Rmu_numubarccoth  = (TH1D*) sk_1Rmu_numubarccoth  -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numubarccoth",  kErecoBinEdges1Rmu);
     h_sk_1Rmu_numubarnc1piC = (TH1D*) sk_1Rmu_numubarnc1piC -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numubarnc1piC", kErecoBinEdges1Rmu);
     h_sk_1Rmu_numubarncoth  = (TH1D*) sk_1Rmu_numubarncoth  -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_numubarncoth",  kErecoBinEdges1Rmu);
     h_sk_1Rmu_nueccqe	     = (TH1D*) sk_1Rmu_nueccqe	     -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nueccqe",       kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuecc1pi      = (TH1D*) sk_1Rmu_nuecc1pi      -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuecc1pi",      kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuecccoh      = (TH1D*) sk_1Rmu_nuecccoh      -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuecccoh",      kErecoBinEdges1Rmu);
     h_sk_1Rmu_nueccoth      = (TH1D*) sk_1Rmu_nueccoth      -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nueccoth",      kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuenc1piC     = (TH1D*) sk_1Rmu_nuenc1piC     -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuenc1piC",     kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuencoth      = (TH1D*) sk_1Rmu_nuencoth      -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuencoth",      kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuebarccqe    = (TH1D*) sk_1Rmu_nuebarccqe    -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuebarccqe",    kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuebarcc1pi   = (TH1D*) sk_1Rmu_nuebarcc1pi   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuebarcc1pi",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuebarcccoh   = (TH1D*) sk_1Rmu_nuebarcccoh   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuebarcccoh",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuebarccoth   = (TH1D*) sk_1Rmu_nuebarccoth   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuebarccoth",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuebarnc1piC  = (TH1D*) sk_1Rmu_nuebarnc1piC  -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuebarnc1piC",  kErecoBinEdges1Rmu);
     h_sk_1Rmu_nuebarncoth   = (TH1D*) sk_1Rmu_nuebarncoth   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_nuebarncoth",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnumubarccqe    = (TH1D*) sk_1Rmu_oscnumubarccqe    -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnumubarccqe",    kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnumubarcc1pi   = (TH1D*) sk_1Rmu_oscnumubarcc1pi   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnumubarcc1pi",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnumubarcccoh   = (TH1D*) sk_1Rmu_oscnumubarcccoh   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnumubarcccoh",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnumubarccoth   = (TH1D*) sk_1Rmu_oscnumubarccoth   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnumubarccoth",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnuebarccqe    = (TH1D*) sk_1Rmu_oscnuebarccqe    -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnuebarccqe",    kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnuebarcc1pi   = (TH1D*) sk_1Rmu_oscnuebarcc1pi   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnuebarcc1pi",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnuebarcccoh   = (TH1D*) sk_1Rmu_oscnuebarcccoh   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnuebarcccoh",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnuebarccoth   = (TH1D*) sk_1Rmu_oscnuebarccoth   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnuebarccoth",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnueccqe    = (TH1D*) sk_1Rmu_oscnueccqe    -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnueccqe",    kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnuecc1pi   = (TH1D*) sk_1Rmu_oscnuecc1pi   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnuecc1pi",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnuecccoh   = (TH1D*) sk_1Rmu_oscnuecccoh   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnuecccoh",   kErecoBinEdges1Rmu);
     h_sk_1Rmu_oscnueccoth   = (TH1D*) sk_1Rmu_oscnueccoth   -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_oscnueccoth",   kErecoBinEdges1Rmu);
  } else {
     h_sk_1Rmu_total         = new TH1D(*sk_1Rmutotal);
     h_sk_1Rmu_numuccqe      = new TH1D(*sk_1Rmu_numuccqe);
     h_sk_1Rmu_numucc1pi     = new TH1D(*sk_1Rmu_numucc1pi);
     h_sk_1Rmu_numucccoh     = new TH1D(*sk_1Rmu_numucccoh);
     h_sk_1Rmu_numuccoth     = new TH1D(*sk_1Rmu_numuccoth);
     h_sk_1Rmu_numunc1piC    = new TH1D(*sk_1Rmu_numunc1piC);
     h_sk_1Rmu_numuncoth     = new TH1D(*sk_1Rmu_numuncoth);
     h_sk_1Rmu_oscnumuccqe      = new TH1D(*sk_1Rmu_oscnumuccqe);
     h_sk_1Rmu_oscnumucc1pi     = new TH1D(*sk_1Rmu_oscnumucc1pi);
     h_sk_1Rmu_oscnumucccoh     = new TH1D(*sk_1Rmu_oscnumucccoh);
     h_sk_1Rmu_oscnumuccoth     = new TH1D(*sk_1Rmu_oscnumuccoth);
     h_sk_1Rmu_numubarccqe   = new TH1D(*sk_1Rmu_numubarccqe);
     h_sk_1Rmu_numubarcc1pi  = new TH1D(*sk_1Rmu_numubarcc1pi);
     h_sk_1Rmu_numubarcccoh  = new TH1D(*sk_1Rmu_numubarcccoh);
     h_sk_1Rmu_numubarccoth  = new TH1D(*sk_1Rmu_numubarccoth);
     h_sk_1Rmu_numubarnc1piC = new TH1D(*sk_1Rmu_numubarnc1piC);
     h_sk_1Rmu_numubarncoth  = new TH1D(*sk_1Rmu_numubarncoth);
     h_sk_1Rmu_nueccqe	     = new TH1D(*sk_1Rmu_nueccqe);	
     h_sk_1Rmu_nuecc1pi      = new TH1D(*sk_1Rmu_nuecc1pi);
     h_sk_1Rmu_nuecccoh      = new TH1D(*sk_1Rmu_nuecccoh);
     h_sk_1Rmu_nueccoth      = new TH1D(*sk_1Rmu_nueccoth);
     h_sk_1Rmu_nuenc1piC     = new TH1D(*sk_1Rmu_nuenc1piC);
     h_sk_1Rmu_nuencoth      = new TH1D(*sk_1Rmu_nuencoth);
     h_sk_1Rmu_nuebarccqe    = new TH1D(*sk_1Rmu_nuebarccqe);
     h_sk_1Rmu_nuebarcc1pi   = new TH1D(*sk_1Rmu_nuebarcc1pi);
     h_sk_1Rmu_nuebarcccoh   = new TH1D(*sk_1Rmu_nuebarcccoh);
     h_sk_1Rmu_nuebarccoth   = new TH1D(*sk_1Rmu_nuebarccoth);
     h_sk_1Rmu_nuebarnc1piC  = new TH1D(*sk_1Rmu_nuebarnc1piC);
     h_sk_1Rmu_nuebarncoth   = new TH1D(*sk_1Rmu_nuebarncoth);
     h_sk_1Rmu_oscnumubarccqe    = new TH1D(*sk_1Rmu_oscnumubarccqe);
     h_sk_1Rmu_oscnumubarcc1pi   = new TH1D(*sk_1Rmu_oscnumubarcc1pi);
     h_sk_1Rmu_oscnumubarcccoh   = new TH1D(*sk_1Rmu_oscnumubarcccoh);
     h_sk_1Rmu_oscnumubarccoth   = new TH1D(*sk_1Rmu_oscnumubarccoth);
     h_sk_1Rmu_oscnueccqe    = new TH1D(*sk_1Rmu_oscnueccqe);
     h_sk_1Rmu_oscnuecc1pi   = new TH1D(*sk_1Rmu_oscnuecc1pi);
     h_sk_1Rmu_oscnuecccoh   = new TH1D(*sk_1Rmu_oscnuecccoh);
     h_sk_1Rmu_oscnueccoth   = new TH1D(*sk_1Rmu_oscnueccoth);
     h_sk_1Rmu_oscnuebarccqe    = new TH1D(*sk_1Rmu_oscnuebarccqe);
     h_sk_1Rmu_oscnuebarcc1pi   = new TH1D(*sk_1Rmu_oscnuebarcc1pi);
     h_sk_1Rmu_oscnuebarcccoh   = new TH1D(*sk_1Rmu_oscnuebarcccoh);
     h_sk_1Rmu_oscnuebarccoth   = new TH1D(*sk_1Rmu_oscnuebarccoth);
  }

//  pdffile->Close();

  //
  // Now re-group pdfs in fewer categories
  //

  TH1D * h_mc_all         = 0;
  TH1D * h_mc_numuccqe    = 0;
  TH1D * h_mc_numuccnonqe = 0;
  TH1D * h_mc_numubarccqe   = 0;
  TH1D * h_mc_numubarccnonqe   = 0;
  TH1D * h_mc_nuenuebarcc = 0;
  TH1D * h_mc_nc          = 0;

  h_mc_all = new TH1D(*h_sk_1Rmu_total);

  h_mc_numuccqe = new TH1D(*h_sk_1Rmu_numuccqe);
  h_mc_numuccqe -> Add (h_sk_1Rmu_oscnumuccqe);

  h_mc_numuccnonqe = new TH1D(*h_sk_1Rmu_numucc1pi);
  h_mc_numuccnonqe -> Add (h_sk_1Rmu_numucccoh);
  h_mc_numuccnonqe -> Add (h_sk_1Rmu_numuccoth);
  h_mc_numuccnonqe -> Add (h_sk_1Rmu_oscnumucc1pi);
  h_mc_numuccnonqe -> Add (h_sk_1Rmu_oscnumucccoh);
  h_mc_numuccnonqe -> Add (h_sk_1Rmu_oscnumuccoth);

  h_mc_numubarccqe   = new TH1D(*h_sk_1Rmu_numubarccqe);
  h_mc_numubarccqe -> Add (h_sk_1Rmu_oscnumubarccqe);

  h_mc_numubarccnonqe  = new TH1D(*h_sk_1Rmu_numubarcc1pi);
  h_mc_numubarccnonqe  -> Add (h_sk_1Rmu_numubarcccoh);
  h_mc_numubarccnonqe  -> Add (h_sk_1Rmu_numubarccoth);
  h_mc_numubarccnonqe  -> Add (h_sk_1Rmu_oscnumubarcc1pi);
  h_mc_numubarccnonqe  -> Add (h_sk_1Rmu_oscnumubarcccoh);
  h_mc_numubarccnonqe  -> Add (h_sk_1Rmu_oscnumubarccoth);

  h_mc_nuenuebarcc  = new TH1D(*h_sk_1Rmu_nueccqe);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_nuecc1pi);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_nuecccoh);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_nueccoth);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_nuebarccqe);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_nuebarcc1pi);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_nuebarcccoh);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_nuebarccoth);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_oscnueccqe);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_oscnuecc1pi);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_oscnuecccoh);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_oscnueccoth);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_oscnuebarccqe);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_oscnuebarcc1pi);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_oscnuebarcccoh);
  h_mc_nuenuebarcc -> Add (h_sk_1Rmu_oscnuebarccoth);
  
  h_mc_nc = new TH1D(*h_sk_1Rmu_numunc1piC);
  h_mc_nc     -> Add (h_sk_1Rmu_numuncoth);
  h_mc_nc     -> Add (h_sk_1Rmu_numubarnc1piC);
  h_mc_nc     -> Add (h_sk_1Rmu_numubarncoth);
  h_mc_nc     -> Add (h_sk_1Rmu_nuenc1piC);
  h_mc_nc     -> Add (h_sk_1Rmu_nuencoth);
  h_mc_nc     -> Add (h_sk_1Rmu_nuebarnc1piC);
  h_mc_nc     -> Add (h_sk_1Rmu_nuebarncoth);

  //
  // If a second pdffile was specified, get the total only
  //

  TH1D * h_mc_all_alt = 0;
  if(altpdffilename) {
    // get pdf
    TFile * altpdffile = new TFile(altpdffilename, "read");
    altpdffile->cd();
    TH1D * sk_1Rmutotal_alt = (TH1D*) pdffile->Get("sk_1Rmutotal");
    // create local copy, rebinned or not
    if(rebin) {
       h_mc_all_alt = (TH1D*) sk_1Rmutotal_alt -> Rebin(kNErecoBins1Rmu, "h_sk_1Rmu_total_alt", kErecoBinEdges1Rmu);
    } else {
       h_mc_all_alt = new TH1D(*sk_1Rmutotal_alt);
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
         h_data = (TH1D*) hd -> Rebin(kNErecoBins1Rmu, "histo_data", kErecoBinEdges1Rmu);
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
  h_mc_numuccnonqe -> SetLineColor(kBlack);
  h_mc_numuccnonqe -> SetFillColor(kWhite);
  h_mc_numuccnonqe -> SetFillStyle(3354);
  h_mc_numuccnonqe -> SetLineColor(kGray+3);
  h_mc_numuccnonqe -> SetFillColor(kGray+3);
  h_mc_numubarccqe   -> SetLineColor(kGreen+1);
  h_mc_numubarccqe   -> SetFillColor(kGreen+1);
  h_mc_numubarccnonqe   -> SetLineColor(kCyan+1);
  h_mc_numubarccnonqe   -> SetFillColor(kCyan+1);
  h_mc_nuenuebarcc -> SetLineColor(kAzure-3);
  h_mc_nuenuebarcc -> SetFillColor(kAzure-3);
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
    h_mc_numuccqe->SetBinContent(ix, (h_mc_numuccqe->GetBinContent(ix) / h_mc_numuccqe->GetBinWidth(ix)) * 0.05);
    h_mc_numuccnonqe->SetBinContent(ix, (h_mc_numuccnonqe->GetBinContent(ix) / h_mc_numuccnonqe->GetBinWidth(ix)) * 0.05);
    h_mc_numubarccqe->SetBinContent(ix, (h_mc_numubarccqe->GetBinContent(ix) / h_mc_numubarccqe->GetBinWidth(ix)) * 0.05);
    h_mc_numubarccnonqe->SetBinContent(ix, (h_mc_numubarccnonqe->GetBinContent(ix) / h_mc_numubarccnonqe->GetBinWidth(ix)) * 0.05);
    h_mc_nuenuebarcc->SetBinContent(ix, (h_mc_nuenuebarcc->GetBinContent(ix) / h_mc_nuenuebarcc->GetBinWidth(ix)) * 0.05);
    if(h_mc_all_alt)
      h_mc_all_alt->SetBinContent(ix, (h_mc_all_alt->GetBinContent(ix) / h_mc_all_alt->GetBinWidth(ix)) * 0.05);
    if(h_data)
      h_data->SetBinContent(ix, (h_data->GetBinContent(ix) / h_data->GetBinWidth(ix)) * 0.05);
  }//ix

  // stack together components of the first pdf
  THStack * stack = new THStack("stack","spectrum contributions");
  stack->Add(h_mc_nc);
  stack->Add(h_mc_nuenuebarcc);
  stack->Add(h_mc_numubarccqe);
  stack->Add(h_mc_numubarccnonqe);
  stack->Add(h_mc_numuccnonqe);
  stack->Add(h_mc_numuccqe);

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
     legend->AddEntry(h_mc_numuccqe,     Form("%s #nu_{#mu} CCQE",pdflabel),            "F");
     legend->AddEntry(h_mc_numuccnonqe,  Form("%s #nu_{#mu} CCnonQE",pdflabel),         "LF");
     legend->AddEntry(h_mc_numubarccqe,    Form("%s #bar{#nu_{#mu}} CCQE",pdflabel),        "LF");
     legend->AddEntry(h_mc_numubarccnonqe,    Form("%s #bar{#nu_{#mu}} CCnonQE",pdflabel),        "LF");
     legend->AddEntry(h_mc_nuenuebarcc,  Form("%s #nu_{e}, #bar{#nu_{e}} CC",pdflabel), "LF");
     legend->AddEntry(h_mc_nc,           Form("%s NC",pdflabel),                        "LF");
  }
  legend->Draw();

  c->SaveAs("spectrum.eps");
  c->SaveAs("spectrum.pdf");
}

