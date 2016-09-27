//____________________________________________________________________________
/*!

  \class     valor::utils

  \brief     Utilities for macros

  \author    Tom Dealtry <t.dealtry \at lancaster.ac.uk>
             Lancaster University
          
  \created   Aug 12, 2015

  \copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _t2xk_macros_utils_h_
#define _t2xk_macros_utils_h_

#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TXMLEngine.h"

#include "valor/Oscillation/OscParams.h"
#include "Enumerations.h"

#include <string>

namespace valor {
  class SuperkExpectation;
  class SuperkObservation;
  class PhysHypothesis;
  class FitMarginalization;
  class RunOpt;

namespace util {
  
  // Tweak the relevant systematics
  void PerformSysThrow(unsigned int tweak_systematics, SuperkExpectation * skexp, bool print = true);
  //TODO should we change the uint to an enum SysTweakId_t?

#ifndef __USE_2D_TEMPLATES__ //1D templates
  void AddExtraEventsToHist(TH1D * h, string * full_points_to_add);
  void AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH1D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand);
#else //2D templates
  void AddExtraEventsToHist(TH2D * h, string * full_points_to_add);
  void AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH2D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand);
#endif
  string * ConstructPointsToAddFromXML(TXMLEngine & xml, XMLNodePointer_t & root_node, valor::conv::SuperkSample_t s, string sextra, int datasetid);

  //void PerformMargFit(int argc, char ** argv, RunOpt *runopt, FitMarginalization *fitmarg, SuperkExpectation *skexp, SuperkObservation *skobs_fakedata, bool print = true); // take inputs from command line
  void PerformMargFit(int argc, char ** argv, RunOpt *runopt, FitMarginalization *fitmarg, SuperkExpectation *skexp, SuperkObservation *skobs_fakedata, bool print = true, string masshier=""); // take inputs from command line

  void MargFit_1D_deltaCP(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_1D_betaBar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_1D_betaBarRateOnly(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_1D_sinsqth13(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_1D_dmsq32(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_1D_dmsq32bar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_2D_sinsqth13_deltaCP(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_2D_deltaCP_sinsqth13(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_2D_sinsqth23_deltaCP(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_2D_dmsq32_deltaCP(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_1D_sinsqth23(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_1D_sinsqth23bar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_2D_sinsqth23_dmsq32(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  //void MargFit_2D_sinsqth23bar_dmsq32bar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_2D_dmsq_dmsqbar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);
  void MargFit_2D_sinsqth23_sinsqth23bar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations);

  // Set the mass hierarchy in Opt and SKMC (for both bar & non-bar) using either
  // - the current value in Opt
  // - a prior
  void SetMassHierarchy(SuperkExpectation * skexp);

  conv::PhysParamId_t GetFreePar(string freeparname);

} // util
} // valor

#endif
