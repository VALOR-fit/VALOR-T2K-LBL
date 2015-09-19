//____________________________________________________________________________
/*!

  \class     valor::utils::OscParams

  \brief     Singleton class to hold the oscillation parameter options

  \author    Tom Dealtry <t.dealtry \at lancaster.ac.uk>
             Lancaster University
          
  \created   Aug 15, 2015

  \copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _valor_osc_params_h_
#define _valor_osc_params_h_

#include "Conventions/NuMassHierarchy.h"
#include "FitCore/PhysHypothesis.h"
#include "FitCore/FitParamList.h"

#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

#include <string>

namespace valor {
namespace util {

class OscParams {
public:

  static OscParams * Instance(void);

  // Check ALL oscillation parameters for -99 (set to global best fits) and -999 (randomize)
  void SetValues(double input_mixing12_nh, double input_mixing13_nh, double input_mixing23_nh,
		 double input_deltacp_nh,  double input_dm12sq_nh,   double input_dm23sq_nh,
		 double input_beta_bar,
		 double input_mixing12_ih, double input_mixing13_ih, double input_mixing23_ih,
		 double input_deltacp_ih,  double input_dm12sq_ih,   double input_dm23sq_ih,
		 double mh_nh_prior,       bool marginalize_MH
		 );
  
  // Randomize ALL oscillation parameters (if requested)
  void Randomize();

  // Update all the values in the PhysHypothesis / FitParamList
  void UpdatePhysHypothesis(PhysHypothesis & ph,  bool print = true);
  void UpdateFitParamList  (FitParamList * fparm, bool set_init, bool print = false);


  // Toy histograms
  void FillOscToyHisto(PhysHypothesis &ph);
  TH1D *GetToyHisto(string oscparam,string MH="normal");
  TH1D *GetToyHisto_Mix13(string MH="normal");
  TH1D *GetToyHisto_Mix23(string MH="normal");
  TH1D *GetToyHisto_Dm23(string MH="normal");
  TH1D *GetToyHisto_Mix12(string MH="normal");
  TH1D *GetToyHisto_Dm12(string MH="normal");
  TH1D *GetToyHisto_dcp(string MH="normal");
  TH1D *GetToyHisto_MH();

  // Write the histograms in the OscParams
  void Write(TFile * f);

  // Print values of oscillation parameters (or whether they're set to be randomized)
  void Print (ostream & stream) const;
  void Print_Toys (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const OscParams & oscparams);

  void SetPrior(conv::PhysParamId_t physpar,string mode,double par1,double par2, string mixmode="");
  
  double ExtractParam(conv::PhysParamId_t physpar,string mode,double par1,double par2,string mixmode="");
  
  string mixmode_mixing12; // "sinsqth" or "sinsq2th"
  string mixmode_mixing13; // "sinsqth" or "sinsq2th"
  string mixmode_mixing23; // "sinsqth" or "sinsq2th"

  string priormode_mixing12; // "user","gaus" or "flat"
  string priormode_mixing13; // "user","gaus" or "flat"
  string priormode_mixing23; // "user","gaus" or "flat"
  string priormode_dcp;      // "user","gaus" or "flat"
  string priormode_dm12sq;   // "user","gaus" or "flat"
  string priormode_dm23sq;   // "user","gaus" or "flat"

  double params_mixing12[2]; // (mean,sigma) or (low,up)
  double params_mixing13[2]; // (mean,sigma) or (low,up)
  double params_mixing23[2]; // (mean,sigma) or (low,up)
  double params_dcp[2];      // (mean,sigma) or (low,up)
  double params_dm12sq[2];   // (mean,sigma) or (low,up) 
  double params_dm23sq[2];   // (mean,sigma) or (low,up)


  //NH
  double mixing12_nh;
  double mixing13_nh;
  double mixing23_nh;
  double dcp_nh;
  double dm12sq_nh;
  double dm23sq_nh;
  //IH
  double mixing12_ih;
  double mixing13_ih;
  double mixing23_ih;
  double dcp_ih;
  double dm12sq_ih;
  double dm23sq_ih;
  //Other
  double beta_bar;
    

  //double GetInitOscParValue(conv::PhysParamId_t oscparid) { return fparm.InitValue(oscparid);}
  //double GetCurrOscParValue(conv::PhysParamId_t oscparid) { return fparm.BestFitValue(oscparid);}

  //Randomize?
  //NH
  bool randomize_mixing12_nh;
  bool randomize_mixing13_nh;
  bool randomize_mixing23_nh;
  bool randomize_dcp_nh;
  bool randomize_dm12sq_nh;
  bool randomize_dm23sq_nh;
  //IH
  bool randomize_mixing12_ih;
  bool randomize_mixing13_ih;
  bool randomize_mixing23_ih;
  bool randomize_dcp_ih;
  bool randomize_dm12sq_ih;
  bool randomize_dm23sq_ih;
  //Other
  bool randomize_beta_bar;

  //Mass hierarchy margninalisation
  double mh_nh_prior;
  int    iMH;
  conv::NuMassHierarchy_t MH;
  bool marginalize_mh;

  //Priors
  std::string priors_file;
  TH2D * hPrior_th13_dcp_nh;
  TH2D * hPrior_th13_dcp_ih;
  TH2D * hPrior_23_nh;
  TH2D * hPrior_23_ih;

  //Save the parameters actually used in histograms
  TH1D *hToyHisto_dcp_nh;
  TH1D *hToyHisto_Mix13_nh;
  TH1D *hToyHisto_Mix23_nh;
  TH1D *hToyHisto_Dm23_nh;
  TH1D *hToyHisto_Mix12_nh;
  TH1D *hToyHisto_Dm12_nh;
  TH1D *hToyHisto_dcp_ih;
  TH1D *hToyHisto_Mix13_ih;
  TH1D *hToyHisto_Mix23_ih;
  TH1D *hToyHisto_Dm23_ih;
  TH1D *hToyHisto_Mix12_ih;
  TH1D *hToyHisto_Dm12_ih;
  TH1D *hToyHisto_MH;

 private:
  // Self
  static OscParams * fInstance;

  void Init();
  void Cleanup();

  // Singleton class: constructors are private
  OscParams();
  OscParams(const OscParams & opt);
  virtual ~OscParams();

  // Clean-up
  struct Cleaner {
    void DummyMethodAndSilentCompiler() { }
    ~Cleaner() {
      if (OscParams::fInstance != 0) {
	delete OscParams::fInstance;
	OscParams::fInstance = 0;
      }
    }
  };
  friend struct Cleaner;
};

ostream & operator << (ostream & stream, const OscParams & oscparams);

} // util
} // valor

#endif
