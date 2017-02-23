//____________________________________________________________________________
/*!
 
 \class    valor::SuperkMC
 
 \brief    SuperK Monte Carlo templates
 
 \author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
 STFC, Rutherford Appleton Laboratory
 
 \created  March 15, 2011
 
 */
//____________________________________________________________________________

#ifndef _valor_sk_mc_templates_h_
#define _valor_sk_mc_templates_h_

#include <TObject.h>
#include <TXMLEngine.h>
#include <TMatrixD.h>
#include <vector>

#include "Conventions.h"
#include "valor/FitCore/PhysHypothesis.h"
#include "valor/Oscillation/NuOsc3Flavour.h"
#include "valor/Oscillation/NuOscSterile.h"
#include "Utilities.h"
#include "valor/VArray/VALORArray.h"


#include <TH3D.h> // NEW 2D TEMPLATES
                  //class TH3D; // NEW 2D TEMPLATES


class TVirtualPad;
class TTree;
class TH1D;
class TH2D;


namespace valor {
  
  class NuFlux;
  class NuXSec;
  class SuperkDetectorFSI;
  
  class SuperkMC: public TObject
  {
    
  public:
    SuperkMC();
    ~SuperkMC();
    
    
    using TObject::Read; // Suppress clang warnings
    using TObject::Draw; // Suppress clang warnings
    
    // Apply systematic variations + oscillations (,decay,decoherence) to MC templates
    void UpdateTemplates(void);

    // Integrate current MC templates and get Ereco spectra (needed for fitting the data)
    bool  UpdateErecoSpectra  (void);
    TH1D* CurrErecoSpectrum   (conv::SuperkSample_t s) const; // Returns new object, mut be deleted
    TH1D* CurrErecoSpectrum   (conv::SuperkSample_t s, conv::SuperkMode_t m) const; // Returns new object, mut be deleted
    double*       CurrErecoSpectrumArray   (conv::SuperkSample_t s);
    const double* CurrErecoSpectrumArray   (conv::SuperkSample_t s) const;
    double*       CurrErecoSpectrumArray   (conv::SuperkSample_t s, conv::SuperkMode_t m);
    const double* CurrErecoSpectrumArray   (conv::SuperkSample_t s, conv::SuperkMode_t m) const;
    VALORArray * CurrSpectrumVArray        (conv::SuperkSample_t s);
    const VALORArray * CurrSpectrumVArray  (conv::SuperkSample_t s) const;
    
    
    
    
    // NEW 2D TEMPLATES
    // Integrate current MC templates and get Ereco Vs Q2 spectra (needed for fitting the data)
    bool  UpdateErecoVsQ2Spectra  (void);
    TH2D* CurrErecoVsQ2Spectrum   (conv::SuperkSample_t s) const; // Returns new object, mut be deleted
    TH2D* CurrErecoVsQ2Spectrum   (conv::SuperkSample_t s, conv::SuperkMode_t m) const; // Returns new object, mut be deleted
    double*       CurrErecoVsQ2SpectrumArray   (conv::SuperkSample_t s);
    const double* CurrErecoVsQ2SpectrumArray   (conv::SuperkSample_t s) const;
    double*       CurrErecoVsQ2SpectrumArray   (conv::SuperkSample_t s, conv::SuperkMode_t m);
    const double* CurrErecoVsQ2SpectrumArray   (conv::SuperkSample_t s, conv::SuperkMode_t m) const;
    VALORArray * CurrSpectrumVArray2D        (conv::SuperkSample_t s);
    const VALORArray * CurrSpectrumVArray2D  (conv::SuperkSample_t s) const;
    //
    
    
    
    // MC event samples used to fill-in MC templates and POT exposure for each
    TTree * MCTree (conv::SuperkSample_t s, conv::SuperkMode_t m);
    double  MCPOT  (conv::SuperkSample_t s, conv::SuperkMode_t m);
    
    // Access MC templates
    TH2D * Template (conv::SuperkSample_t s, conv::SuperkMode_t m, conv::SKMCTemplateType_t type) const;
    TH2D * VarTemplate (conv::SuperkSample_t s, conv::SuperkMode_t m) const;
    
    TH2D Get2DHist(conv::SuperkSample_t s, conv::SKMCTemplateType_t type) const;
    TH2D Get2DHist(conv::SuperkSample_t s, conv::SuperkMode_t m , conv::SKMCTemplateType_t type) const;
    
    
    
    // NEW 2D TEMPLATES
    TH3D * Template2D (conv::SuperkSample_t s, conv::SuperkMode_t m, conv::SKMCTemplateType_t type) const;
    TH3D * VarTemplate2D (conv::SuperkSample_t s, conv::SuperkMode_t m) const;
    TH3D Get3DHist(conv::SuperkSample_t s, conv::SKMCTemplateType_t type) const;
    TH3D Get3DHist(conv::SuperkSample_t s, conv::SuperkMode_t m , conv::SKMCTemplateType_t type) const;
    //
    
    
    
    // POT of analysed data sample (MC is scaled to that POT value)
    double AnalyzedSamplePOT(conv::SuperkSample_t s) const;
    
    // Read inputs from XML file
    bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
    bool ReadSmearedMomentum(TXMLEngine* xml, XMLNodePointer_t root_node);
    
    // Physics hypothesis
    PhysHypothesis & Hypothesis(void);
    
    // Hooks for obtaining necessary flux, cross-section and detector information
    void TellMeAboutTheFlux       (NuFlux *         flux )       { fNuFlux = flux; }
    void TellMeAboutCrossSections (NuXSec *         xsec )       { fNuXSec = xsec; }
    void TellMeAboutTheDetectorFSI(SuperkDetectorFSI * skdetfsi) { fSKDetFSI  = skdetfsi; }
    
    // Calculate error in each 2D template bin due to finite MC statistics & fluctuate to study effect of MC stat err
    void CalculateMCStatError(void); // called at init, using the original templates
    void ApplyMCStatFluctationToTemplates(bool shape_only=true); // fluctuate contents of MC templates using MC stat error
    
    // Draw inputs:
    // ...
    void Draw(TVirtualPad * c, Option_t * opt = "") const;
    
    // Hack for x_mec until we have a parameter we can tweak
    void SetX2p2h(double xmec)              { fX2p2h = xmec; }
    
    // Allow osc prob > 1 in 2 flavour oscillation?
    void SetUnphysicalOscProb(bool unphys) { fUnphysicalOscProb = unphys; }
    
    // Allow changing the mass hierarchy
    void SetMassHierarchy       (conv::NuMassHierarchy_t hier) { fNuOsc3Flavour.SetMassHierarchy(hier); }
    // NEW
    void SetMassHierarchyBar    (conv::NuMassHierarchy_t hier) { fNuOsc3FlavourBar.SetMassHierarchy(hier); }
    void SetMassHierarchySterile(conv::NuMassHierarchy_t hier) { fNuOscSterile.SetMassHierarchy(hier); }
    //
    
    
    void NevDebug(conv::SKMCTemplateType_t t, string name);
    void NevDebug2D(conv::SKMCTemplateType_t t, string name); // NEW 2D TEMPLATES
    
    void FillTemplateArrays(bool with_tuning_weights = false);
    bool DoTemplateArraysHaveNIWGTune(void) const;
    
    // NEW 2D TEMPLATES
    
    double CalcThetaDegrees(float dir[10][3]);
    
    double CalcQ2(double Enu,double Elep,double Plep,double costheta,double Mlep);
    double CalcQd2(double Enu,double Elep,double Plep,double costheta,double Mlep);
    double CalcEhad(double Enu,double Elep);
    double CalcInelasticity(double Enu,double Elep);
    double CalcBjorkenScale(double Enu,double Elep,double Plep,double costheta,double Mlep,double Mn);
    double CalcWhadr(double Enu,double Elep,double Plep,double costheta,double Mlep,double Mn);
    double CalcXSecTerm(double Enu,double Elep,double Plep,double costheta,double Mlep,double Mn);
    double ErecoFunc_Nu(double Mom,double costh,double mass_lep);
    double ErecoFunc_NuBar(double Mom,double costh,double mass_lep);
    //
    
    double SmearMomentum(double fit_momentum,Option_t * opt = "") const;
    
    
  private:
    
    // Init and clean-up
    void Init    (void);
    void CleanUp (void);
    void Book    (void);
    
    // Fill MC templates from the input SK event tree or from other sources
    void FillTemplateHistogram(conv::SuperkSample_t s, conv::SuperkMode_t m);
    void FillTemplateHistogramSmeared(conv::SuperkSample_t s, conv::SuperkMode_t m);
    void FillSpecialTemplates (void);

    
    // NEW 2D TEMPLATES
    void DrawProjectionXY(TH3D *h);  // for debugging
    void ProjectTemplate3Dto2D(TH3D *h3d[conv::kns][conv::knm], TH2D *h2d[conv::kns][conv::knm]);
    double GetTemplateIntegral(TH3D *h3d[conv::kns][conv::knm]);
    double GetTemplateIntegral(TH2D *h2d[conv::kns][conv::knm]);
    double GetTemplateIntegral(TH3D *h3d);
    double GetTemplateIntegral(TH2D *h2d);
    //
    
    // Apply BANFF flux + cross section tuning (based on ND280 data) on top of NA61-based flux tuned MC
    void ApplyBANFFTuningToTemplates (void);
    void ApplyBANFFTuningToTemplates2D (void); // NEW 2D TEMPLATES
    
    // Apply systematic variations to MC templates
    void ApplySystematicsToTemplates                   (void);
    void ApplySystematicsToTemplates_Multiplicative    (void);
    void ApplySystematicsToTemplates_NonMultiplicative (void);
    bool ApplyEnergyScaleSystematicToTotalSpectrum     (void);
    
    // NEW 2D TEMPLATES
    void ApplySystematicsToTemplates2D                 (void);
    void ApplySystematicsToTemplates_Multiplicative_2D (void);
    void ApplySystematicsToTemplates_NonMultiplicative2D (void);
    bool ApplyEnergyScaleSystematicToTotalSpectrum2D   (void);
    //
    
    // Apply physics hypothesis (oscillations,...) to MC templates
    void ApplyPhysHypothesisToTemplates (void); // Master function, picks one of the following functions to call based on fPhysHypothesis
    
    // NEW 2D TEMPLATES
    void ApplyPhysHypothesisToTemplates2D (void); // Master function, picks one of the following functions to call based on fPhysHypothesis
                                                  //
    
    
    // NEW 2D TEMPLATES
    /* void ApplyPhysHypothesisToTemplates_3fOsc(void);  */
    /* void ApplyPhysHypothesisToTemplates_3fOscBar(void);  */
    /* void ApplyPhysHypothesisToTemplates_3fOscBarNormBar(void);   */
    /* void ApplyPhysHypothesisToTemplates_3fOscNormBar(void);      */
    /* void ApplyPhysHypothesisToTemplates_Sterile(void);  */
    void ApplyPhysHypothesisToTemplates_3fOsc(const bool Use2Dtemplates=false);
    void ApplyPhysHypothesisToTemplates_3fOscBar(const bool Use2Dtemplates=false);
    void ApplyPhysHypothesisToTemplates_3fOscBarDissap(const bool Use2Dtemplates=false);
    
    void ApplyPhysHypothesisToTemplates_3fOscBarNormBar(const bool Use2Dtemplates=false);
    void ApplyPhysHypothesisToTemplates_3fOscNormBar(const bool Use2Dtemplates=false);
    void ApplyPhysHypothesisToTemplates_Sterile(const bool Use2Dtemplates=false);
    //
    
    
    // Ignore nue_signal templates for all but the 3-flavour oscillation hypothesis
    bool SkipTemplate(conv::SuperkMode_t m);
    
    // Debug method, prints out some number of events, and also saves templates & spectra
    // void NevDebug(conv::SKMCTemplateType_t t, string name);
    void MCStatErrDebug(void);
    
    // Looks to see if there are any events for the given Etrue in the sample
    bool IsThereAnEvent(conv::SuperkSample_t s, double Etrue) const;
    
    // Integrate current MC templates and get Etrue spectra (needed for applying physics hypothesis)
    void FillEtrueSpectra(void);
    void FillEtrueSpectra2D(void); // NEW 2D TEMPLATES
    
    //
    // Private data members
    //
    
    double fPOT;           ///< POT for actual data-set (MC is scaled to that exposure)
    double fPOTRHC;        ///< POT for actual data-set in RHC sample (MC is scaled to that exposure)
    double fPOTMultiRing;  ///< POT for actual data-set in FHC MultiRing sample (MC is scaled to that exposure)
    
    // input MC templates (nominal, nominal + systematics, nominal + systematics + physics hypothesis) and Ereco spectra
    TH2D *     fMCTemplOri[conv::kns][conv::knm];        ///< Ereco-Etrue MC template(sample,mode): nominal (doesn't even have pot weighting applied)
    TH2D *     fMCTemplNom[conv::kns][conv::knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning
    TH2D *     fMCTemplRPA[conv::kns][conv::knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning + with rel RPA
    TH2D *     fMCTemplPOT[conv::kns][conv::knm];        ///< Ereco-Etrue MC template(sample,mode): nominal + ext/tuning + with rel RPA
    
    TH2D *     fMCVariableTempl[conv::kns][conv::knm];
    
    std::vector<VALORArray * > fMCVariableTemplArray; // sample -> (mode, MCetrue,Xsereco)
    std::vector<VALORArray * > fMCTemplArrayOri;      // sample -> (mode, etrue, ereco)
    
    std::vector<VALORArray * > fMCTemplArrayNom; //sample -> (mode, etrue, observable)
    std::vector<VALORArray * > fMCTemplArrayTun; //sample -> (mode, etrue, observable)
    std::vector<VALORArray * > fMCTemplArrayFin; //sample -> (mode, etrue, observable)
    
    std::vector<VALORArray * > fMCTemplArrayRPAWght; //sample -> (mode, etrue, observable) - contains weight RPA_rel / weight RPA_nonrel
    
    std::vector<VALORArray * > fMCStatFracErr;           // sample -> (mode, etrue, observable)
    std::vector<VALORArray * > fMCEtrueSpectrumTotArray; // sample -> etrue
    std::vector<VALORArray * > fMCErecoSpectrumTotArray; // sample -> observable
    std::vector<VALORArray * > fMCErecoSpectrumArray   ; // sample -> (mode, observable)
    
    TH1D *     fMCErecoSpectrum[conv::kns][conv::knm];   ///< Ereco MC spectrum (sample,mode)     : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
    TH1D *     fMCErecoSpectrumTot[conv::kns];     ///< Ereco MC spectrum (sample)          : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
    TH1D *     fMCEtrueSpectrumTot[conv::kns];     ///< Etrue MC spectrum (sample)          : nominal
    
    
    
    
    
    // NEW 2D TEMPLATES
    TH3D *     fMCTemplOri2D[conv::kns][conv::knm];        ///< Ereco-Q2-Etrue MC template(sample,mode): nominal (doesn't even have pot weighting applied)
    TH3D *     fMCTemplNom2D[conv::kns][conv::knm];        ///< Ereco-Q2-Etrue MC template(sample,mode): nominal + ext/tuning
    TH3D *     fMCTemplRPA2D[conv::kns][conv::knm];        ///< Ereco-Q2-Etrue MC template(sample,mode): nominal + ext/tuning + with rel RPA
    TH3D *     fMCTemplPOT2D[conv::kns][conv::knm];        ///< Ereco-Q2-Etrue MC template(sample,mode): nominal + ext/tuning + with rel RPA
    
    TH3D *     fMCVariableTempl2D[conv::kns][conv::knm];
    
    std::vector<VALORArray * > fMCVariableTemplArray2D; // sample -> (mode, MCetrue,Xsereco,XsQ2)
    std::vector<VALORArray * > fMCTemplArrayOri2D;      // sample -> (mode, etrue, ereco, Q2)
    
    std::vector<VALORArray * > fMCTemplArrayNom2D; //sample -> (mode, etrue, observable1, observable2)
    std::vector<VALORArray * > fMCTemplArrayTun2D; //sample -> (mode, etrue, observable1, observable2)
    std::vector<VALORArray * > fMCTemplArrayFin2D; //sample -> (mode, etrue, observable1, observable2)
    
    std::vector<VALORArray * > fMCTemplArrayRPAWght2D; //sample -> (mode, etrue, observable1, observable2) - contains weight RPA_rel / weight RPA_nonrel
    
    std::vector<VALORArray * > fMCStatFracErr2D;               // sample -> (mode, etrue, observable1, observable2)
                                                               //std::vector<VALORArray * > fMCEtrueSpectrumTotArray2D;   // sample -> etrue
    std::vector<VALORArray * > fMCErecoSpectrumTotArray2D;     // sample -> observable1
    std::vector<VALORArray * > fMCQ2SpectrumTotArray2D;        // sample -> observable2
    std::vector<VALORArray * > fMCErecoVsQ2SpectrumTotArray2D; // sample -> observable1, observable2
    std::vector<VALORArray * > fMCErecoVsQ2SpectrumArray;      // sample -> (mode, observable1, observable2)
    
    TH2D *     fMCErecoVsQ2Spectrum[conv::kns][conv::knm];   ///< Ereco Vs Q2 MC spectrum (sample,mode)     : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
    TH2D *     fMCErecoVsQ2SpectrumTot[conv::kns];     ///< Ereco Vs Q2 MC spectrum (sample)          : nominal + ext/tuning + nd280/tuning + systematics + physics hypothesis
                                                       //
    
    
    // input event trees and normalisation factors
    vector<TTree*> fMCTreeNuMu;              ///< input MC event tree /  sample: numu                            / cuts: by sample
    vector<TTree*> fMCTreeNuMuBar;           ///< input MC event tree /  sample: numubar                         / cuts: by sample
    vector<TTree*> fMCTreeNuE;               ///< input MC event tree /  sample: nue (intrinsic)                 / cuts: by sample
    vector<TTree*> fMCTreeNuEBar;            ///< input MC event tree /  sample: nuebar (intrinsic)              / cuts: by sample
    vector<TTree*> fMCTreeOscNuE;            ///< input MC event tree /  sample: nue (signal, 100% numu->nue)    / cuts: by sample
    vector<TTree*> fMCTreeOscNuEBar;         ///< input MC event tree /  sample: nuebar (signal, 100% numu->nue) / cuts: by sample
    
    // FHC
    TH1D*      fMCPOTHistNuMu;               ///< input MC event sample POT histogram : numu
    TH1D*      fMCPOTHistNuMuBar;            ///< input MC event sample POT histogram : numubar
    TH1D*      fMCPOTHistNuE;                ///< input MC event sample POT histogram : nue
    TH1D*      fMCPOTHistNuEBar;             ///< input MC event sample POT histogram : nuebar
    TH1D*      fMCPOTHistOscNuE;             ///< input MC event sample POT histogram : nue (signal)
    TH1D*      fMCPOTHistOscNuEBar;          ///< input MC event sample POT histogram : nuebar (signal)
    double     fMCPOTNuMu;                   ///< MC event sample POT : numu
    double     fMCPOTNuMuBar;                ///< MC event sample POT : numubar
    double     fMCPOTNuE;                    ///< MC event sample POT : nue
    double     fMCPOTNuEBar;                 ///< MC event sample POT : nuebar
    double     fMCPOTOscNuE;                 ///< MC event sample POT : nue (signal)
    double     fMCPOTOscNuEBar;              ///< MC event sample POT : nuebar (signal)
    
    // RHC
    TH1D*      fMCPOTHistRHCNuMu;            ///< input MC event sample POT histogram : numu
    TH1D*      fMCPOTHistRHCNuMuBar;         ///< input MC event sample POT histogram : numubar
    TH1D*      fMCPOTHistRHCNuE;             ///< input MC event sample POT histogram : nue
    TH1D*      fMCPOTHistRHCNuEBar;          ///< input MC event sample POT histogram : nuebar
    TH1D*      fMCPOTHistRHCOscNuE;          ///< input MC event sample POT histogram : nue (signal)
    TH1D*      fMCPOTHistRHCOscNuEBar;       ///< input MC event sample POT histogram : nuebar (signal)
    double     fMCPOTRHCNuMu;                ///< MC event sample POT : numu
    double     fMCPOTRHCNuMuBar;             ///< MC event sample POT : numubar
    double     fMCPOTRHCNuE;                 ///< MC event sample POT : nue
    double     fMCPOTRHCNuEBar;              ///< MC event sample POT : nuebar
    double     fMCPOTRHCOscNuE;              ///< MC event sample POT : nue (signal)
    double     fMCPOTRHCOscNuEBar;           ///< MC event sample POT : nuebar (signal)
    
    // MultiRing FHC
    TH1D*      fMCPOTHistMultiRingNuMu;               ///< input MC event sample POT histogram : numu
    TH1D*      fMCPOTHistMultiRingNuMuBar;            ///< input MC event sample POT histogram : numubar
    TH1D*      fMCPOTHistMultiRingNuE;                ///< input MC event sample POT histogram : nue
    TH1D*      fMCPOTHistMultiRingNuEBar;             ///< input MC event sample POT histogram : nuebar
    TH1D*      fMCPOTHistMultiRingOscNuE;             ///< input MC event sample POT histogram : nue (signal)
    TH1D*      fMCPOTHistMultiRingOscNuEBar;          ///< input MC event sample POT histogram : nue (signal)
    double     fMCPOTMultiRingNuMu;                   ///< MC event sample POT : numu
    double     fMCPOTMultiRingNuMuBar;                ///< MC event sample POT : numubar
    double     fMCPOTMultiRingNuE;                    ///< MC event sample POT : nue
    double     fMCPOTMultiRingNuEBar;                 ///< MC event sample POT : nuebar
    double     fMCPOTMultiRingOscNuE;                 ///< MC event sample POT : nue (signal)
    double     fMCPOTMultiRingOscNuEBar;              ///< MC event sample POT : nuebar (signal)
    
    ////
    float      fBrPtrue[10];                 ///< branch var for reading input MC Minituple
    double     fBrEreco;                     ///< branch var for reading input MC mini-DST
                                             //~ double     fBrEtrue;                     ///< branch var for reading input MC mini-DST
    int        fBrGenCode;                   ///< branch var for reading input MC mini-DST
                                             //~ double     fBrWght;                      ///< branch var for reading input MC mini-DST
                                             //~ float      fBrWght1Re;                   ///< branch var for reading input MC mini-DST
    double     fBrWghtFlx;                   ///< branch var for reading input MC mini-DST
    double     fBrWghtSFRFG;                 ///< branch var for reading input MC mini-DST
    double     fBrWghtRPArel;                ///< branch var for reading input MC mini-DST
    double     fBrWghtRPAnonrel;             ///< branch var for reading input MC mini-DST
    double     fBrWghtCCcoh;                 ///< branch var for reading input MC mini-DST
    double     fBrWghtOsc;                   ///< branch var for reading input MC mini-DST, used for comparison with p-theta analysis oscillation weights as branch
    int        fBrEvFlag;                    ///< branch var for reading input MC mini-DST
                                             // NEW 2D TEMPLATES
                                             //int        fBrNumnu;                     ///< branch var for reading input MC mini-DST
                                             //int        fBrNring;                     ///< branch var for reading input MC mini-DST
    float      fBrAmomm[10];                  ///< branch var for reading input MC mini-DST (1 mu-like ring)
    float      fBrAmome[10];                  ///< branch var for reading input MC mini-DST (1 e-like ring)
    float      fBrDir[1][3];                    ///< branch var for reading input MC mini-DST
    float      fBrDirNu[10][3];                 ///< branch var for reading input MC mini-DST
                                                //
    
    static const int kMaxVarBins   = 1000;
    // NEW 2D TEMPLATES
    static const int kMaxVarBins2D = 10000;
    //
    
    // Arrays for remembering which big bins have events, and which little bins for those big bigs have events.
    int fNNonEmptyVariableBins[conv::kns][conv::knm];               //!
    int fNNonEmptyFixedBins   [conv::kns][conv::knm][kMaxVarBins];  //!
    int fNonEmptyVariableBins [conv::kns][conv::knm][kMaxVarBins];  //! These contain both the reco and true bin numbers, masked with conv::kBigIndex;
    int * fNonEmptyFixedBins  [conv::kns][conv::knm][kMaxVarBins];  //! These contain both the reco and true bin numbers, masked with conv::kBigIndex;
    
    
    // NEW 2D TEMPLATES
    int fNNonEmptyVariableBins2D[conv::kns][conv::knm];               //!
    int fNNonEmptyFixedBins2D   [conv::kns][conv::knm][kMaxVarBins2D];  //!
    int fNonEmptyVariableBins2D [conv::kns][conv::knm][kMaxVarBins2D];  //! These contain both the q2, reco and true bin numbers, masked with conv::kBigIndex;
    int * fNonEmptyFixedBins2D  [conv::kns][conv::knm][kMaxVarBins2D];  //! These contain both the q2, reco and true bin numbers, masked with conv::kBigIndex;
                                                                        //
    
    
    // utility objects for computing cross-section and flux weights
    NuFlux *         fNuFlux;                ///< Neutrino flux inputs
    NuXSec *         fNuXSec;                ///< Neutrino cross-section inputs
    SuperkDetectorFSI * fSKDetFSI;           ///< SuperK detector efficiency error etc
    
    // utility objects for defining physics hypothesis (oscillations, decay, ...) & computing probabilities
    PhysHypothesis fPhysHypothesis;          ///< Physics hypothesis
    NuOsc3Flavour  fNuOsc3Flavour;           ///< 3-flavour neutrino osc. probability calculations
    NuOsc3Flavour  fNuOsc3FlavourBar;        ///< 3-flavour anti-neutrino osc. probability calculations
    
    // NEW
    NuOscSterile  fNuOscSterile;            ///< sterile osc. probability calculations
                                            //
    
    double fX2p2h; //hack x_mec until we have a parameter we can use
    
    bool fUnphysicalOscProb; // Allow unphysical region (probs not in range [0,1]) in 2 flavour case
    
    bool fTemplateArraysHaveNIWGTune;
    
    util::ArrIdx * arr_index; //! Pointer to singleton, do not save in TObject
    
    ClassDef(SuperkMC, 8);
  };
  
} //valor

#endif
