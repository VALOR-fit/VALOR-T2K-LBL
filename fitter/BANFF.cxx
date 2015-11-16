//________________________________________________________________________________________
/*

 \copyright (c) 2010-2013, VaLOR group.

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__  
//#define __show_init_mesg_debug__   
//#define __show_cholesky_mesg_info__
//#define __show_cholesky_mesg_debug__
//#define __show_tweak_mesg_debug__
//#define __show_update_correlated_mesg_debug__
//#define __show_chisqpenalty_mesg_debug__
//#define __show_chisqpenalty_mesg_debug2__
//#define __BANFF_BESTFIT_IS_ONE__ //debug mode to turn the best fit vector into a vector of 1's (i.e. BANFF tuned == flux tuned)
//#define __show_banffsystconv_check_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip> 
#include <TRandom3.h>
#include <TMath.h>

#include "Build.h"
#include "Enumerations.h"
#include "Utilities.h"
#include "BANFF.h" 
#include "Conventions.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace valor;
using namespace valor::conv;
using namespace valor::util;

ClassImp(BANFF)

//________________________________________________________________________________________
BANFF::BANFF() :
TObject()
{
  this->Init();
}
//________________________________________________________________________________________
BANFF::~BANFF()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void BANFF::Init(void)
{
  fBestFitValues        = 0; 
  fCurrentValues        = 0; 
  fCovMtrx              = 0;
  fCovMtrxInv           = 0;
  fCovMtrxChoL          = 0; 

  fRHCBestFitValues     = 0;
  fRHCCurrentValues     = 0;
  fRHCCovMtrx           = 0;
  fRHCCovMtrxInv        = 0;
  fRHCCovMtrxChoL       = 0;
  
  // NEW FLUX FHC+RHC
  fFHCpRHCBestFitValues     = 0;
  fFHCpRHCCurrentValues     = 0;
  fFHCpRHCCovMtrx           = 0;
  fFHCpRHCCovMtrxInv        = 0;
  fFHCpRHCCovMtrxChoL       = 0;
  //

  fIsTweaked            = false;
  fIsFluxWghtUpdated    = true;
  fIsTweakedRHC         = false;
  fIsFluxWghtUpdatedRHC = true;
  // NEW FLUX FHC+RHC
  //fIsTweakedFHCpRHC         = false;
  //fIsFluxWghtUpdatedFHCpRHC = true;
  //
  
  // NEW FHC+RHC JOINT  
  // for (int ip = 0 ; ip < conv::knppmax ; ip++) {
  //   PhysParamId_t p = (PhysParamId_t) ip;
  //   fParam2Idx[ip]=-1;
  //   for (int ibp = 0 ; ibp < kNumOfParamsInBANFFfit ; ibp++) {
  //     if (kListOfBANFFParams[ibp] == p) {
  //       fParam2Idx[ip] = ibp;
  //     }
  //   }
  //   for (int ibp = 0 ; ibp < kNumOfRHCParamsInBANFFfit ; ibp++) {
  //     if (kListOfRHCBANFFParams[ibp] == p) {
  //       fParam2Idx[ip] = ibp;
  //     }
  //   }
  // }

  
  for (int ip = 0 ; ip < conv::knppmax ; ip++) {
    PhysParamId_t p = (PhysParamId_t) ip;
    // FHC-only
    fParam2Idx_FHC[ip]=-1;
    for (int ibp = 0 ; ibp < kNumOfParamsInBANFFfit ; ibp++) {
      if (kListOfBANFFParams[ibp] == p) {
        fParam2Idx_FHC[ip] = ibp;
      }
    }
    // RHC-only
    fParam2Idx_RHC[ip]=-1;
    for (int ibp = 0 ; ibp < kNumOfRHCParamsInBANFFfit ; ibp++) {
      if (kListOfRHCBANFFParams[ibp] == p) {
        fParam2Idx_RHC[ip] = ibp;
      }
    }
    // FHC+RHC 
    fParam2Idx_FHCpRHC[ip]=-1;
    for (int ibp = 0 ; ibp < kNumOfFHCpRHCParamsInBANFFfit ; ibp++) {
      if (kListOfFHCpRHCBANFFParams[ibp] == p) {
        fParam2Idx_FHCpRHC[ip] = ibp;
      }
    }
  }
  //
  
  fUntunedCentralValues = 0;
  fRHCUntunedCentralValues = 0;
  fFHCpRHCUntunedCentralValues = 0;   // NEW FHC+RHC JOINT

}
//________________________________________________________________________________________
void BANFF::CleanUp(void)
{
  if ( fBestFitValues             ) delete fBestFitValues;
  if ( fCurrentValues             ) delete fCurrentValues;
  if ( fCovMtrx                   ) delete fCovMtrx;
  if ( fCovMtrxChoL               ) delete fCovMtrxChoL;
  
  if ( fRHCBestFitValues          ) delete fRHCBestFitValues;
  if ( fRHCCurrentValues          ) delete fRHCCurrentValues;
  if ( fRHCCovMtrx                ) delete fRHCCovMtrx;
  if ( fRHCCovMtrxChoL            ) delete fRHCCovMtrxChoL;

  // NEW FLUX FHC+RHC
  if ( fFHCpRHCBestFitValues          ) delete fFHCpRHCBestFitValues;
  if ( fFHCpRHCCurrentValues          ) delete fFHCpRHCCurrentValues;
  if ( fFHCpRHCCovMtrx                ) delete fFHCpRHCCovMtrx;
  if ( fFHCpRHCCovMtrxChoL            ) delete fFHCpRHCCovMtrxChoL;
  //
  
  if (fUntunedCentralValues)        delete fUntunedCentralValues;
  if (fRHCUntunedCentralValues)     delete fRHCUntunedCentralValues;
  if (fFHCpRHCUntunedCentralValues) delete fFHCpRHCUntunedCentralValues; // NEW FHC+RHC JOINT
}
//________________________________________________________________________________________
bool BANFF::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read BANFF inputs to oscillation analysis
//
  if (fBestFitValues)           {delete fBestFitValues;           fBestFitValues           = 0 ; }
  if (fCurrentValues)           {delete fCurrentValues;           fCurrentValues           = 0 ; }
  if (fUntunedCentralValues)    {delete fUntunedCentralValues;    fUntunedCentralValues    = 0 ; } 
  //if (fRHCUntunedCentralValues) {delete fRHCUntunedCentralValues; fRHCUntunedCentralValues = 0 ; } // NEW FHC+RHC JOINT
  if (fCovMtrx)                 {delete fCovMtrx;                 fCovMtrx                 = 0 ; }
  if (fCovMtrxInv)              {delete fCovMtrxInv;              fCovMtrxInv              = 0 ; }
  
  const char * banff_param_path = T2KOpt::Instance()->fPathBanffParam.c_str();
#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Reading BANFF inputs from XML configuration file." << endl;
  cout << "Exp/BANFF >> "<<util::XML2RootFileName(xml,root_node,banff_param_path)
       << " - " << util::XML2RootObjName (xml,root_node,banff_param_path);
#endif

  const int np = conv::kNumOfParamsInBANFFfit;
  
#ifdef __USE_NIWG2015__
  
  // NEW FHC+RHC JOINT
  fUntunedCentralValues    = util::XML2TVectorD(xml,root_node,"xsec/untuned_parameter_central_values");
  if (fUntunedCentralValues->GetNoElements()!=kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit) {
    cerr << "NuXSec/Read >> Expected "<<kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit
         << " elements in <untuned_parameter_central_values>, got "<<fUntunedCentralValues->GetNoElements()
         << endl;
    exit(-1);
  }
  //

  // fUntunedCentralValues    = util::XML2TVectorD(xml,root_node,"xsec/untuned_parameter_central_values");
  // fRHCUntunedCentralValues = util::XML2TVectorD(xml,root_node,"xsec/untuned_parameter_central_values"); 
  // if (not (fUntunedCentralValues && fRHCUntunedCentralValues) ) {
  // cerr<<"Failed to read <untuned_parameter_central_values>"<<endl;
  // exit(-1);
  // }
  // if (fUntunedCentralValues->GetNoElements()!=kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit) {
  // cerr << "NuXSec/Read >> Expected "<<kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit
  // << " elements in <untuned_parameter_central_values>, got "<<fUntunedCentralValues->GetNoElements()
  // << endl;
  // exit(-1);
  // }
  // if (fRHCUntunedCentralValues->GetNoElements()!=kNumOfRHCParamsInBANFFfit-kNumOfRHCFluxParamsInBANFFfit) {
  // cerr << "NuXSec/Read >> Expected "<<kNumOfRHCParamsInBANFFfit-kNumOfRHCFluxParamsInBANFFfit
  // << " elements in <untuned_parameter_central_values>, got "<<fRHCUntunedCentralValues->GetNoElements()
  // << endl;
  // exit(-1);
  // }
  
#endif



  // Get BANFF fit central values

#ifndef __BANFF_BESTFIT_IS_ONE__
  fBestFitValues = util::XML2TVectorD(xml,root_node,banff_param_path);
  if(!fBestFitValues) {
    cerr << "Exp/BANFF >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
#else
  double init[np];
  for(int ip = 0; ip < np; ip++)
    init[ip] = 1;
  fBestFitValues = new TVectorD(np, init);
#endif

  const TVectorD & bestfit = *fBestFitValues;
  if(bestfit.GetNoElements() != np) 
  {
    cerr << "Exp/BANFF >> BANFF best-fit parameter vector has unexpected number of elements!" << endl;
    cerr << "Exp/BANFF >> Expecting: " <<  np << " elements" << endl;
    cerr << "Exp/BANFF >> Found: " << bestfit.GetNoElements() << " elements" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> Input BANFF best-fit parameters: " << endl; 
  util::PrintCompact(bestfit,12,6,"BANFF fit params");
#endif

  // init current values to best-fit values
  fCurrentValues = new TVectorD(bestfit);

  // Get BANFF fit covariance matrix
  fCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathBanffCovMatr.c_str());
 
  if(!fCovMtrx) {
    cerr << "Exp/BANFF >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  // Invert BANFF fit covariance matrix
  double det;
  fCovMtrxInv = new TMatrixD(*fCovMtrx);
  (*fCovMtrxInv).Invert(&det);

  const TMatrixD & cov = *fCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np) 
  {
    cerr << "Exp/BANFF >> BANFF covariance matrix has unexpected size!" << endl;
    cerr << "Exp/BANFF >> Expecting: " <<  np << " x " << np << " matrix" << endl;
    cerr << "Exp/BANFF >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
  }
#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> Input BANFF covariance matrix: " << endl; 
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF >> Inverse BANFF covariance matrix: " << endl;
  util::PrintCompact(*fCovMtrxInv,7,5);
#endif

#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Done reading BANFF inputs!" << endl;
#endif
  
  fBestFitArray = fBestFitValues->GetMatrixArray();
  fCurrentArray = fCurrentValues->GetMatrixArray();
  return true;
}
//________________________________________________________________________________________
bool BANFF::ReadRHC(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  if (fRHCBestFitValues)         { delete fRHCBestFitValues; fRHCBestFitValues = 0; }
  if (fRHCCurrentValues)         { delete fRHCCurrentValues; fRHCCurrentValues = 0; }
  if (fRHCCovMtrx)               { delete fRHCCovMtrx;       fRHCCovMtrx       = 0; }
  if (fRHCCovMtrxInv)            { delete fRHCCovMtrxInv;    fRHCCovMtrxInv    = 0; }
  // NEW FHC+RHC JOINT
  if (fRHCUntunedCentralValues)   {delete fRHCUntunedCentralValues; fRHCUntunedCentralValues = 0 ; }
  //

  // Read BANFF inputs to oscillation analysis
  //
#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Reading RHC BANFF inputs from XML configuration file." << endl;
#endif

  const int np = conv::kNumOfRHCParamsInBANFFfit;


  // NEW FHC+RHC JOINT
#ifdef __USE_NIWG2015__
  fRHCUntunedCentralValues = util::XML2TVectorD(xml,root_node,"xsec/untuned_parameter_central_values");
  if (fRHCUntunedCentralValues->GetNoElements()!=kNumOfRHCParamsInBANFFfit-kNumOfRHCFluxParamsInBANFFfit) {
    cerr << "NuXSec/Read >> Expected "<<kNumOfRHCParamsInBANFFfit-kNumOfRHCFluxParamsInBANFFfit
	 << " elements in <untuned_parameter_central_values>, got "<<fRHCUntunedCentralValues->GetNoElements()
	 << endl;
    exit(-1);
  }
#endif
  //


  // Get BANFF fit central values

#ifndef __BANFF_BESTFIT_IS_ONE__

  fRHCBestFitValues = util::XML2TVectorD(xml,root_node,T2KOpt::Instance()->fPathRHCBanffParam.c_str());
    
  if(!fRHCBestFitValues) {
    cerr << "Exp/BANFF >> No RHC BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
#else
  double init[np];
  for(int ip = 0; ip < np; ip++)
    init[ip] = 1;
  fRHCBestFitValues = new TVectorD(np, init);
#endif

  const TVectorD & bestfit = *fRHCBestFitValues;
  if(bestfit.GetNoElements() != np) 
    {
      cerr << "Exp/BANFF >> RHC BANFF best-fit parameter vector has unexpected number of elements!" << endl;
      cerr << "Exp/BANFF >> Expecting: " <<  np << " elements" << endl;
      cerr << "Exp/BANFF >> Found: " << bestfit.GetNoElements() << " elements" << endl;
      exit(1);
    }

#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> RHC input BANFF best-fit parameters: " << endl; 
  util::PrintCompact(bestfit,12,6,"RHC BANFF fit params");
#endif

  // init current values to best-fit values
  fRHCCurrentValues = new TVectorD(bestfit);

  // Get BANFF fit covariance matrix
  
  fRHCCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathRHCBanffCovMatr.c_str());

  if(!fRHCCovMtrx) {
    cerr << "Exp/BANFF >> No RHC BANFF covariance matrix was found." << endl;
    exit(1);
  }
  // Invert BANFF fit covariance matrix
  double det;
  fRHCCovMtrxInv = new TMatrixD(*fRHCCovMtrx);
  (*fRHCCovMtrxInv).Invert(&det);

  const TMatrixD & cov = *fRHCCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np) 
    {
      cerr << "Exp/BANFF >> RHC BANFF covariance matrix has unexpected size!" << endl;
      cerr << "Exp/BANFF >> Expecting: " <<  np << " x " << np << " matrix" << endl;
      cerr << "Exp/BANFF >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
    }
#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> RHC input BANFF covariance matrix: " << endl; 
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF >> RHC inverse BANFF covariance matrix: " << endl;
  util::PrintCompact(*fRHCCovMtrxInv,7,5);
#endif

#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Done reading RHC BANFF inputs!" << endl;
#endif
  
  fRHCBestFitArray = fRHCBestFitValues->GetMatrixArray();
  fRHCCurrentArray = fRHCCurrentValues->GetMatrixArray();
  return true;
}
//________________________________________________________________________________________



// NEW FLUX FHC+RHC

bool BANFF::ReadFHCpRHC(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  if (fBestFitValues)  { delete fBestFitValues; fBestFitValues = 0;}
  if (fCurrentValues)  { delete fCurrentValues; fCurrentValues = 0;}
  if (fCovMtrx)        { delete fCovMtrx;       fCovMtrx       = 0;}
  if (fCovMtrxInv)     { delete fCovMtrxInv;    fCovMtrxInv    = 0;}

  if (fRHCBestFitValues)  { delete fRHCBestFitValues; fRHCBestFitValues = 0;}
  if (fRHCCurrentValues)  { delete fRHCCurrentValues; fRHCCurrentValues = 0;}
  if (fRHCCovMtrx)        { delete fRHCCovMtrx;       fRHCCovMtrx       = 0;}
  if (fRHCCovMtrxInv)     { delete fRHCCovMtrxInv;    fRHCCovMtrxInv    = 0;}

  if (fFHCpRHCBestFitValues)  { delete fFHCpRHCBestFitValues; fFHCpRHCBestFitValues = 0;}
  if (fFHCpRHCCurrentValues)  { delete fFHCpRHCCurrentValues; fFHCpRHCCurrentValues = 0;}
  if (fFHCpRHCCovMtrx)        { delete fFHCpRHCCovMtrx;       fFHCpRHCCovMtrx       = 0;}
  if (fFHCpRHCCovMtrxInv)     { delete fFHCpRHCCovMtrxInv;    fFHCpRHCCovMtrxInv    = 0;}

  // NEW FHC+RHC JOINT
  if (fUntunedCentralValues)        {delete fUntunedCentralValues;        fUntunedCentralValues    = 0 ; }
  if (fRHCUntunedCentralValues)     {delete fRHCUntunedCentralValues;     fRHCUntunedCentralValues = 0 ; }
  if (fFHCpRHCUntunedCentralValues) {delete fFHCpRHCUntunedCentralValues; fFHCpRHCUntunedCentralValues = 0 ; }
  //
  
  // Read BANFF inputs to oscillation analysis
  //
#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Reading FHC+RHC BANFF inputs from XML configuration file." << endl;
#endif
  const int np = conv::kNumOfFHCpRHCParamsInBANFFfit;



  
  // NEW FHC+RHC JOINT 

  // f*UntunedCentralValues contain only xsec parameters,the same for FHC-only, RHC-only and FHC+RHC
  // All these vectors have the same dimension
#ifdef __USE_NIWG2015__
  fFHCpRHCUntunedCentralValues    = util::XML2TVectorD(xml,root_node,"xsec/untuned_parameter_central_values");
  if (fFHCpRHCUntunedCentralValues->GetNoElements()!=kNumOfFHCpRHCParamsInBANFFfit-kNumOfFHCpRHCFluxParamsInBANFFfit) {
    cerr << "NuXSec/Read >> Expected "<<kNumOfFHCpRHCParamsInBANFFfit-kNumOfFHCpRHCFluxParamsInBANFFfit
	 << " elements in <untuned_parameter_central_values>, got "<<fFHCpRHCUntunedCentralValues->GetNoElements()
	 << endl;
    exit(-1);
  }
  fUntunedCentralValues    = util::XML2TVectorD(xml,root_node,"xsec/untuned_parameter_central_values");
  if (fUntunedCentralValues->GetNoElements()!=kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit) {
    cerr << "NuXSec/Read >> Expected "<<kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit
	 << " elements in <untuned_parameter_central_values>, got "<<fUntunedCentralValues->GetNoElements()
	 << endl;
    exit(-1);
  }
  fRHCUntunedCentralValues    = util::XML2TVectorD(xml,root_node,"xsec/untuned_parameter_central_values");
  if (fRHCUntunedCentralValues->GetNoElements()!=kNumOfRHCParamsInBANFFfit-kNumOfRHCFluxParamsInBANFFfit) {
    cerr << "NuXSec/Read >> Expected "<<kNumOfRHCParamsInBANFFfit-kNumOfRHCFluxParamsInBANFFfit
	 << " elements in <untuned_parameter_central_values>, got "<<fRHCUntunedCentralValues->GetNoElements()
	 << endl;
    exit(-1);
  }
#endif
  //



  // Get BANFF fit central values

#ifndef __BANFF_BESTFIT_IS_ONE__

  // Use same path as FHC. Also FHC and RHC paths are the same. Use different configuration files   
  fFHCpRHCBestFitValues = util::XML2TVectorD(xml,root_node,T2KOpt::Instance()->fPathBanffParam.c_str());
  
  cout << T2KOpt::Instance()->fPathBanffParam.c_str() << endl; 
  
  if(!fFHCpRHCBestFitValues) {
    cerr << "Exp/BANFF >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
#else
  double init[np];
  for(int ip = 0; ip < np; ip++)
    init[ip] = 1;
  fFHCpRHCBestFitValues = new TVectorD(np, init);
#endif

  const TVectorD & bestfit = *fFHCpRHCBestFitValues;
  if(bestfit.GetNoElements() != np) 
  {
    cerr << "Exp/BANFF >> FHC+RHC BANFF best-fit parameter vector has unexpected number of elements!" << endl;
    cerr << "Exp/BANFF >> Expecting: " <<  np << " elements" << endl;
    cerr << "Exp/BANFF >> Found: " << bestfit.GetNoElements() << " elements" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> Input FHC+RHC BANFF best-fit parameters: " << endl; 
  util::PrintCompact(bestfit,12,6,"FHC+RHC BANFF fit params");
#endif  
  
  // init current values to best-fit values
  fFHCpRHCCurrentValues = new TVectorD(bestfit);

  // Get BANFF fit covariance matrix

  // Use same path as FHC. Also FHC and RHC paths are the same. Use different configuration files 
  fFHCpRHCCovMtrx = util::XML2TMatrixD(xml,root_node,T2KOpt::Instance()->fPathBanffCovMatr.c_str());
 
  if(!fFHCpRHCCovMtrx) {
    cerr << "Exp/BANFF >> No FHC+RHC BANFF covariance matrix was found." << endl;
    exit(1);
  }
  // Invert BANFF fit covariance matrix
  double det;
  fFHCpRHCCovMtrxInv = new TMatrixD(*fFHCpRHCCovMtrx);
  (*fFHCpRHCCovMtrxInv).Invert(&det);


  // TODO: get the FHC-only and RHC-only covariance matrices!!!
  // write a method to select the sub-matrices from the full FHC+RHC one
  fCovMtrx = new TMatrixD(conv::kNumOfParamsInBANFFfit,conv::kNumOfParamsInBANFFfit);
  fCovMtrxInv = new TMatrixD(conv::kNumOfParamsInBANFFfit,conv::kNumOfParamsInBANFFfit);
  fRHCCovMtrx = new TMatrixD(conv::kNumOfRHCParamsInBANFFfit,conv::kNumOfRHCParamsInBANFFfit);
  fRHCCovMtrxInv = new TMatrixD(conv::kNumOfRHCParamsInBANFFfit,conv::kNumOfRHCParamsInBANFFfit);
  fCovMtrx->UnitMatrix();
  fCovMtrxInv->UnitMatrix();
  fRHCCovMtrx->UnitMatrix();
  fRHCCovMtrxInv->UnitMatrix();
  //


  const TMatrixD & cov = *fFHCpRHCCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np) 
  {
    cerr << "Exp/BANFF >> FHC+RHC BANFF covariance matrix has unexpected size!" << endl;
    cerr << "Exp/BANFF >> Expecting: " <<  np << " x " << np << " matrix" << endl;
    cerr << "Exp/BANFF >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
  }
#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> Input FHC+RHC BANFF covariance matrix: " << endl; 
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF >> Inverse FHC+RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(*fFHCpRHCCovMtrxInv,7,5);
#endif

#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Done reading BANFF inputs!" << endl;
#endif
  
  fFHCpRHCBestFitArray = fFHCpRHCBestFitValues->GetMatrixArray();
  fFHCpRHCCurrentArray = fFHCpRHCCurrentValues->GetMatrixArray();
  

  // NEW FHC+RHC JOINT

  // Initialize FHC-only and RHC-only vectors

  const int np_fhc = conv::kNumOfParamsInBANFFfit;
  double init_fhc[np_fhc];
  for(int ip = 0; ip < np_fhc; ip++) init_fhc[ip] = conv::kBadResult;
  fBestFitValues = new TVectorD(np_fhc, init_fhc);
  TVectorD & bestfit_values_fhc = *fBestFitValues;

  const int np_rhc = conv::kNumOfRHCParamsInBANFFfit;
  double init_rhc[np_rhc];
  for(int ip = 0; ip < np_rhc; ip++) init_rhc[ip] = conv::kBadResult;
  fRHCBestFitValues = new TVectorD(np_rhc, init_rhc);
  TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  
  this->SplitFHCpRHCvectors(bestfit,bestfit_values_fhc,bestfit_values_rhc);
  
  // Copy bestfit to current values
  fCurrentValues = new TVectorD(bestfit_values_fhc);
  fRHCCurrentValues = new TVectorD(bestfit_values_rhc);

  fBestFitArray = fBestFitValues->GetMatrixArray();
  fCurrentArray = fCurrentValues->GetMatrixArray();
  fRHCBestFitArray = fRHCBestFitValues->GetMatrixArray();
  fRHCCurrentArray = fRHCCurrentValues->GetMatrixArray();

#ifdef __show_init_mesg_debug__
  cout << "Exp/BANFF >> Input FHC BANFF best-fit parameters from FHC+RHC: " << endl; 
  util::PrintCompact(bestfit_values_fhc,12,6,"FHC BANFF fit params");
  cout << "Exp/BANFF >> Input RHC BANFF best-fit parameters from FHC+RHC: " << endl; 
  util::PrintCompact(bestfit_values_rhc,12,6,"RHC BANFF fit params");
#endif

#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Done reading BANFF inputs!" << endl;
#endif
  //
 
  return true;
}

//

//________________________________________________________________________________________
double BANFF::BestFitValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx_FHC(p);
  return fBestFitArray[idx];
}
//________________________________________________________________________________________
double BANFF::BestFitValue(int idx) const
{
  return fBestFitArray[idx];
}
//________________________________________________________________________________________
double BANFF::RHCBestFitValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx_RHC(p);
  return fRHCBestFitArray[idx];
}
//________________________________________________________________________________________
double BANFF::RHCBestFitValue(int idx) const
{
  return fRHCBestFitArray[idx];
}
//________________________________________________________________________________________


// NEW FLUX FHC+RHC

// Not really useful. For Debugging
double BANFF::FHCpRHCBestFitValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx_FHCpRHC(p);
  return fFHCpRHCBestFitArray[idx];
}
//________________________________________________________________________________________
double BANFF::FHCpRHCBestFitValue(int idx) const
{
  return fFHCpRHCBestFitArray[idx];
}
//


//________________________________________________________________________________________
double BANFF::CurrentValue(int idx) const
{
  return fCurrentArray[idx];
}
//________________________________________________________________________________________
double BANFF::CurrentValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx_FHC(p); // NEW FHC+RHC JOINT
  //int idx = this->Param2Idx(p); 
  return fCurrentArray[idx];
}
//________________________________________________________________________________________
double BANFF::RHCCurrentValue(int idx) const
{
  return fRHCCurrentArray[idx];
}
//________________________________________________________________________________________
double BANFF::RHCCurrentValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx_RHC(p); // NEW FHC+RHC JOINT
  //int idx = this->Param2Idx(p); 
  return fRHCCurrentArray[idx];
}
//________________________________________________________________________________________


// NEW FLUX FHC+RHC

// Not really useful. For Debugging
double BANFF::FHCpRHCCurrentValue(int idx) const
{
  return fFHCpRHCCurrentArray[idx];
}
//________________________________________________________________________________________
double BANFF::FHCpRHCCurrentValue(PhysParamId_t p) const
{
  int idx = this->Param2Idx_FHCpRHC(p);
  return fFHCpRHCCurrentArray[idx];
}
//________________________________________________________________________________________
//


double BANFF::OneSigmaError(PhysParamId_t p) const
{
  const TMatrixD & cov = *fCovMtrx;
  int idx = this->Param2Idx_FHC(p); // NEW FHC+RHC JOINT
  //int idx = this->Param2Idx(p);
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double BANFF::OneSigmaError(int idx) const
{
  const TMatrixD & cov = *fCovMtrx;
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double BANFF::RHCOneSigmaError(PhysParamId_t p) const
{
  const TMatrixD & cov = *fRHCCovMtrx;
  int idx = this->Param2Idx_RHC(p); // NEW FHC+RHC JOINT
  //int idx = this->Param2Idx(p); 
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double BANFF::RHCOneSigmaError(int idx) const
{
  const TMatrixD & cov = *fRHCCovMtrx;
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________


// NEW FLUX FHC+RHC
double BANFF::FHCpRHCOneSigmaError(PhysParamId_t p) const
{
  const TMatrixD & cov = *fFHCpRHCCovMtrx;
  int idx = this->Param2Idx_FHCpRHC(p);
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double BANFF::FHCpRHCOneSigmaError(int idx) const
{
  const TMatrixD & cov = *fFHCpRHCCovMtrx;
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________


// NEW FHC+RHC JOINT

//________________________________________________________________________________________
double BANFF::SampleOneSigmaError(PhysParamId_t p) const
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) 
  return this->OneSigmaError(p);
  #else
  return this->FHCpRHCOneSigmaError(p);
  #endif
  //
#elif defined(__FHC_RUNNING__)
  return this->OneSigmaError(p);
#elif defined(__RHC_RUNNING__)
  return this->RHCOneSigmaError(p);
#else
  std::cerr<<"BANFF/SampleOneSigmaError >> Cannot find right matrix."<<endl;
  exit(-1);
#endif
}
//________________________________________________________________________________________
double BANFF::SampleOneSigmaError(int idx) const
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) 
  return this->OneSigmaError(idx);
  #else
  return this->FHCpRHCOneSigmaError(idx);
  #endif
  //
#elif defined(__FHC_RUNNING__)
  return this->OneSigmaError(idx);
#elif defined(__RHC_RUNNING__)
  return this->RHCOneSigmaError(idx);
#else
  std::cerr<<"BANFF/SampleOneSigmaError >> Cannot find right matrix."<<endl;
  exit(-1);
#endif
}
//________________________________________________________________________________________
double BANFF::SampleCurrentValue(PhysParamId_t p) const
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) 
  return this->CurrentValue(p);
  #else
  return this->FHCpRHCCurrentValue(p);
  #endif
  //
#elif defined(__FHC_RUNNING__)
  return this->CurrentValue(p);
#elif defined(__RHC_RUNNING__)
  return this->RHCCurrentValue(p);
#else
  std::cerr<<"BANFF/SampleCurrentValue >> Cannot find right matrix."<<endl;
  exit(-1);
#endif
}
//________________________________________________________________________________________
double BANFF::SampleCurrentValue(int idx) const
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) 
  return this->CurrentValue(idx);
  #else
  return this->FHCpRHCCurrentValue(idx);
  #endif
  //
#elif defined(__FHC_RUNNING__)
  return this->CurrentValue(idx);
#elif defined(__RHC_RUNNING__)
  return this->RHCCurrentValue(idx);
#else
  std::cerr<<"BANFF/SampleCurrentValue >> Cannot find right matrix."<<endl;
  exit(-1);
#endif
}
//________________________________________________________________________________________
double BANFF::SampleBestFitValue(PhysParamId_t p) const
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) 
  return this->BestFitValue(p);
  #else
  return this->FHCpRHCBestFitValue(p);
  #endif
  //
#elif defined(__FHC_RUNNING__)
  return this->BestFitValue(p);
#elif defined(__RHC_RUNNING__)
  return this->RHCBestFitValue(p);
#else
  std::cerr<<"BANFF/SampleBestFitValue >> Cannot find right matrix."<<endl;
  exit(-1);
#endif
}
//________________________________________________________________________________________
double BANFF::SampleBestFitValue(int idx) const
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) 
  return this->BestFitValue(idx);
  #else
  return this->FHCpRHCBestFitValue(idx);
  #endif
  //
#elif defined(__FHC_RUNNING__)
  return this->BestFitValue(idx);
#elif defined(__RHC_RUNNING__)
  return this->RHCBestFitValue(idx);
#else
  std::cerr<<"BANFF/SampleBestFitValue >> Cannot find right matrix."<<endl;
  exit(-1);
#endif
}
//________________________________________________________________________________________

//

//________________________________________________________________________________________
void BANFF::SyncRHCToFHCSystematics()
{
  this->TweakRHC(true, fCurrentValues);
}
//________________________________________________________________________________________


void BANFF::TweakSample(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
  #if defined(__FULLY_CORRELATED__)
  this->Tweak(force_inputs, input_vec);
  this->SyncRHCToFHCSystematics();
  #else
  this->TweakFHCpRHC(force_inputs, input_vec);
  #endif
  //
#else
#ifdef __FHC_RUNNING__
  this->Tweak(force_inputs, input_vec);
#endif
#ifdef __RHC_RUNNING__
  this->TweakRHC(force_inputs, input_vec);
#endif
#endif
}



//________________________________________________________________________________________
void BANFF::Tweak(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking BANFF fit parameters" << endl;
#endif

  // Switch flags
  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //
  if(fCovMtrxChoL == 0) 
    this->CholeskyDecomposeMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fBestFitValues;
  TVectorD & current_values = *fCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
  {
     // Input vector is NULL.
     // This resets the BANFF params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset BANFF params to best-fit values." << endl;
#endif
        fIsTweaked = false;
        current_values = bestfit_values;
        return;
     }
     else     
     {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Forcing input BANFF parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;
        assert(vec.GetNoElements() == conv::kNumOfParamsInBANFFfit);
        current_values = vec;
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/Tweak >> Generating random twk dials (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
     current_values = vec;
     return;
  }

}
//________________________________________________________________________________________
void BANFF::TweakRHC(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking RHC BANFF fit parameters" << endl;
#endif

  // Switch flags
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdatedRHC = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fRHCBestFitValues;
  TVectorD & current_values = *fRHCCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset RHC BANFF params to best-fit values." << endl;
#endif
        fIsTweakedRHC = false;
        current_values = bestfit_values;
        return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/Tweak >> Forcing RHC input BANFF parameters. " << endl;
#endif
	  TVectorD & vec = *input_vec;
	  assert(vec.GetNoElements() == conv::kNumOfRHCParamsInBANFFfit);
	  current_values = vec;
	  return;
	}
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/Tweak >> Generating RHC random twk dials (correlated)! " << endl;
#endif
      TMatrixD & L = *fRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
      current_values = vec;
      return;
    }

}
//________________________________________________________________________________________




// NEW FLUX FHC+RHC
void BANFF::TweakFHCpRHC(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking FHC+RHC BANFF fit parameters" << endl;
#endif

  // Switch flags
  
  // NEW FHC+RHC JOINT
  fIsTweaked            = true;
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdated    = false;
  fIsFluxWghtUpdatedRHC = false;
  //fIsTweakedFHCpRHC         = true;
  //fIsFluxWghtUpdatedFHCpRHC = false;
  //

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //
  
  if(fFHCpRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeFHCpRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  // NEW FHC+RHC JOINT
  TVectorD & bestfit_values_fhc = *fBestFitValues;
  TVectorD & current_values_fhc = *fCurrentValues;
  TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  TVectorD & current_values_rhc = *fRHCCurrentValues;
  //

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset FHC+RHC BANFF params to best-fit values." << endl;
#endif

	// NEW FHC+RHC JOINT
        fIsTweaked    = false;
        fIsTweakedRHC = false;
	//fIsTweakedFHCpRHC = false;
	
        current_values = bestfit_values;
	// NEW FHC+RHC JOINT
        current_values_fhc = bestfit_values_fhc;
        current_values_rhc = bestfit_values_rhc;
        //

	// bestfit_values.Print();
	// bestfit_values_fhc.Print();
	// bestfit_values_rhc.Print();

	// current_values.Print();
	// current_values_fhc.Print();
	// current_values_rhc.Print();
	
	return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/Tweak >> Forcing FHC+RHC input BANFF parameters. " << endl;
#endif
	  TVectorD & vec = *input_vec;
	  assert(vec.GetNoElements() == conv::kNumOfFHCpRHCParamsInBANFFfit);
	  current_values = vec;

	  // NEW FHC+RHC JOINT
	  this->SplitFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);
	  //
	  
	  return;
	}
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/Tweak >> Generating FHC+RHC random twk dials (correlated)! " << endl;
#endif

      TMatrixD & L = *fFHCpRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
      current_values = vec;

      // NEW FHC+RHC JOINT
      this->SplitFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);
      //

      return;
    }
}

//________________________________________________________________________________________

void BANFF::TweakSampleConstrainedOnly(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec   , /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
  bool   unconstrained     /* This function will either do all constrained parameters, or all unconstrained*/
)
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
  #if defined(__FULLY_CORRELATED__)
  this->TweakConstrainedOnly(force_inputs, input_vec,unconstrained);
  this->SyncRHCToFHCSystematics();
  #else
  this->TweakFHCpRHCConstrainedOnly(force_inputs, input_vec,unconstrained);
  #endif
  //
#else
#ifdef __FHC_RUNNING__
  this->TweakConstrainedOnly(force_inputs, input_vec,unconstrained);
#endif
#ifdef __RHC_RUNNING__
  this->TweakRHCConstrainedOnly(force_inputs, input_vec,unconstrained);
#endif
#endif
}



//________________________________________________________________________________________
void BANFF::TweakConstrainedOnly(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec   ,/* input vector of BANFF params - Set to 0 to reset to best-fit ones */
  bool       unconstrained
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking BANFF fit parameters" << endl;
#endif

  // Switch flags
  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //
  if(fCovMtrxChoL == 0) 
    this->CholeskyDecomposeMatrix();

  //
  // Get best-fit and current BANFF params
  //

  
  TVectorD & bestfit_values = *fBestFitValues;
  TVectorD & current_values = *fCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
  {
     // Input vector is NULL.
     // This resets the BANFF params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset BANFF params to best-fit values." << endl;
#endif
        fIsTweaked = false;
        current_values = bestfit_values;
        return;
     }
     else     
     {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Forcing input BANFF parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;
        assert(vec.GetNoElements() == conv::kNumOfParamsInBANFFfit);
        for (int ip = 0; ip < vec.GetNoElements() ; ip++) {
          PhysParamId_t p = conv::kListOfBANFFParams[ip];
          if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec[ip]=current_values[ip];};
        }
        current_values = vec;
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/Tweak >> Generating random twk dials (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
      for (int ip = 0; ip < vec.GetNoElements() ; ip++) {
        PhysParamId_t p = conv::kListOfBANFFParams[ip];
        if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec[ip]=current_values[ip];};
      }
     current_values = vec;
     return;
  }

}
//________________________________________________________________________________________
void BANFF::TweakRHCConstrainedOnly(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec    ,/* input vector of BANFF params - Set to 0 to reset to best-fit ones */
 bool       unconstrained
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking RHC BANFF fit parameters" << endl;
#endif

  // Switch flags
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdatedRHC = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fRHCBestFitValues;
  TVectorD & current_values = *fRHCCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset RHC BANFF params to best-fit values." << endl;
#endif
        fIsTweakedRHC = false;
        current_values = bestfit_values;
        return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/Tweak >> Forcing RHC input BANFF parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;
        assert(vec.GetNoElements() == conv::kNumOfRHCParamsInBANFFfit);
        for (int ip = 0; ip < vec.GetNoElements() ; ip++) {
          PhysParamId_t p = conv::kListOfRHCBANFFParams[ip];
          if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec[ip]=current_values[ip];};
        }
        current_values = vec;
        return;
      }
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/Tweak >> Generating RHC random twk dials (correlated)! " << endl;
#endif
      TMatrixD & L = *fRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
        for (int ip = 0; ip < vec.GetNoElements() ; ip++) {
        PhysParamId_t p = conv::kListOfRHCBANFFParams[ip];
        if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec[ip]=current_values[ip];};
      }
      
      current_values = vec;
      return;
    }

}
//________________________________________________________________________________________




// NEW FLUX FHC+RHC
void BANFF::TweakFHCpRHCConstrainedOnly(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec   , /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
 bool       unconstrained
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking FHC+RHC BANFF fit parameters" << endl;
#endif

  // Switch flags

  // NEW FHC+RHC JOINT
  fIsTweaked            = true;
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdated    = false;
  fIsFluxWghtUpdatedRHC = false;
  //fIsTweakedFHCpRHC         = true;
  //fIsFluxWghtUpdatedFHCpRHC = false;
  //

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fFHCpRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeFHCpRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  // NEW FHC+RHC JOINT
  TVectorD & bestfit_values_fhc = *fBestFitValues;
  TVectorD & current_values_fhc = *fCurrentValues;
  TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  TVectorD & current_values_rhc = *fRHCCurrentValues;
  //

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset FHC+RHC BANFF params to best-fit values." << endl;
#endif
        // NEW FHC+RHC JOINT
	fIsTweaked    = false;
        fIsTweakedRHC = false;
        //fIsTweakedFHCpRHC = false;
	//
	
        current_values = bestfit_values;
	// NEW FHC+RHC JOINT
        current_values_fhc = bestfit_values_fhc;
        current_values_rhc = bestfit_values_rhc;
        //
        return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/Tweak >> Forcing FHC+RHC input BANFF parameters. " << endl;
#endif

	  
	  // 
	  // OLD
	  //
	  // TVectorD & vec = *input_vec;
	  // assert(vec.GetNoElements() == conv::kNumOfFHCpRHCParamsInBANFFfit);
	  // for (int ip = 0; ip < vec.GetNoElements() ; ip++) {
	  //   PhysParamId_t p = conv::kListOfFHCpRHCBANFFParams[ip];
	  //   if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec[ip]=current_values[ip];};
	  // }
	  // current_values = vec;
	  // // NEW FHC+RHC JOINT
	  // this->SplitFHCpRHCvectors(fFHCpRHCBestFitValues,fBestFitValues,fRHCBestFitValues);
	  // //


	  TVectorD & vec = *input_vec;
	  assert(vec.GetNoElements() == conv::kNumOfFHCpRHCParamsInBANFFfit);
	  // Split into fhc-only and rhc-only
	  const int np_fhc = conv::kNumOfParamsInBANFFfit;
	  const int np_rhc = conv::kNumOfRHCParamsInBANFFfit;     
	  TVectorD vec_fhc(np_fhc);
	  TVectorD vec_rhc(np_rhc);
	  this->SplitFHCpRHCvectors(vec,vec_fhc,vec_rhc);      
	  // Copy to fhc-only vectors      
	  for (int ip = 0; ip < vec_fhc.GetNoElements() ; ip++) {
	    PhysParamId_t p = conv::kListOfBANFFParams[ip];
	    if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec_fhc[ip]=current_values_fhc[ip];};
	  }
	  current_values_fhc = vec_fhc;
	  // Copy to rhc-only vectors      
	  for (int ip = 0; ip < vec_rhc.GetNoElements() ; ip++) {
	    PhysParamId_t p = conv::kListOfRHCBANFFParams[ip];
	    if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec_rhc[ip]=current_values_rhc[ip];};
	  }
	  current_values_rhc = vec_rhc;

	  return;
	}
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/Tweak >> Generating RHC random twk dials (correlated)! " << endl;
#endif

      //
      // OLD
      //
      //TMatrixD & L = *fFHCpRHCCovMtrxChoL;
      //TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);    
      // Copy to fhc-only and rhc-only vectors
      // for (int ip = 0; ip < vec.GetNoElements() ; ip++) {
      // 	PhysParamId_t p = conv::kListOfFHCpRHCBANFFParams[ip];
      //   if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec[ip]=current_values[ip];};
      // }
      // current_values = vec;
      // // NEW FHC+RHC JOINT
      // this->SplitFHCpRHCvectors(fFHCpRHCCurrentValues,fCurrentValues,fRHCCurrentValues);
      // //

      TMatrixD & L = *fFHCpRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);    
      // Split into fhc-only and rhc-only
      const int np_fhc = conv::kNumOfParamsInBANFFfit;
      const int np_rhc = conv::kNumOfRHCParamsInBANFFfit;     
      TVectorD vec_fhc(np_fhc);
      TVectorD vec_rhc(np_rhc);
      this->SplitFHCpRHCvectors(vec,vec_fhc,vec_rhc);      
      // Copy to fhc-only vectors      
      for (int ip = 0; ip < vec_fhc.GetNoElements() ; ip++) {
	PhysParamId_t p = conv::kListOfBANFFParams[ip];
	if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec_fhc[ip]=current_values_fhc[ip];};
      }
      current_values_fhc = vec_fhc;
      // Copy to rhc-only vectors      
      for (int ip = 0; ip < vec_rhc.GetNoElements() ; ip++) {
	PhysParamId_t p = conv::kListOfRHCBANFFParams[ip];
	if (unconstrained!=util::IsBANFFUnconstrained(p)) {vec_rhc[ip]=current_values_rhc[ip];};
      }
      current_values_rhc = vec_rhc;

      // vec.Print();
      // fCurrentValues->Print();
      // fRHCCurrentValues->Print();
      // exit(1);

      return;
    }
}


//________________________________________________________________________________________
void BANFF::TweakSampleXSecOnly(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
  #if defined(__FULLY_CORRELATED__)
  this->TweakXSecOnly(force_inputs,input_vec);
  this->SyncRHCToFHCSystematics();
  #else
  this->TweakFHCpRHCXSecOnly(force_inputs,input_vec);
  #endif
  //
#else
#ifdef __FHC_RUNNING__
  this->TweakXSecOnly(force_inputs,input_vec);
#endif
#ifdef __RHC_RUNNING__
  this->TweakRHCXSecOnly(force_inputs,input_vec);
#endif
#endif
}

//________________________________________________________________________________________
void BANFF::TweakXSecOnly(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/TweakXSecOnly >> Tweaking BANFF xsec parameters" << endl;
#endif

  // Switch flags
  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the banff covariance matrix
  //

  if(fCovMtrxChoL == 0) 
    this->CholeskyDecomposeMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fBestFitValues;
  TVectorD & current_values = *fCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
  {
     // Input vector is NULL.
     // This resets the BANFF params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakXSecOnly >> Reset BANFF xsec params to best-fit values." << endl;
#endif
        fIsTweaked = false;
        current_values = bestfit_values;
        return;
     }
     else     
     {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakXSecOnly >> Forcing input BANFF xsec parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;

	int knumofxsecparamsinbanfffit = conv::kNumOfParamsInBANFFfit - conv::kNumOfFluxParamsInBANFFfit; 
        assert(vec.GetNoElements() == knumofxsecparamsinbanfffit);
	for (int i=0; i< knumofxsecparamsinbanfffit; i++){
	  current_values[i] = vec[i];
	}
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/TweakFluxOnly >> Generating random twk dials (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);

     for (int i=0; i< conv::kNumOfParamsInBANFFfit; i++){
       if(i>=conv::kNumOfFluxParamsInBANFFfit) current_values[i] = vec[i]; // xsec params are after flux params
       else current_values[i] = bestfit_values[i];
     }
     return;
  }
}
//________________________________________________________________________________________
void BANFF::TweakRHCXSecOnly(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/TweakXSecOnly >> Tweaking RHC BANFF xsec parameters" << endl;
#endif

  // Switch flags
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdatedRHC = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the banff covariance matrix
  //

  if(fRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fRHCBestFitValues;
  TVectorD & current_values = *fRHCCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
  {
     // Input vector is NULL.
     // This resets the BANFF params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakXSecOnly >> Reset RHC BANFF xsec params to best-fit values." << endl;
#endif
        fIsTweakedRHC = false;
        current_values = bestfit_values;
        return;
     }
     else     
     {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakXSecOnly >> Forcing input RHC BANFF xsec parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;

	int knumofxsecparamsinbanfffit = conv::kNumOfRHCParamsInBANFFfit - conv::kNumOfRHCFluxParamsInBANFFfit; 
        assert(vec.GetNoElements() == knumofxsecparamsinbanfffit);
	for (int i=0; i< knumofxsecparamsinbanfffit; i++){
	  current_values[i] = vec[i];
	}
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/TweakFluxOnly >> Generating RHC random twk dials (correlated)! " << endl;
#endif
     TMatrixD & L = *fRHCCovMtrxChoL;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);

     for (int i=0; i< conv::kNumOfRHCParamsInBANFFfit; i++){
       if(i>=conv::kNumOfRHCFluxParamsInBANFFfit) current_values[i] = vec[i]; // xsec params are after flux params
       else current_values[i] = bestfit_values[i];
     }
     return;
  }
}
//________________________________________________________________________________________
void BANFF::TweakFHCpRHCXSecOnly(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/TweakXSecOnly >> Tweaking FHCpRHC BANFF xsec parameters" << endl;
#endif

  // Switch flags

  // NEW FHC+RHC JOINT
  fIsTweaked            = true;
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdated    = false;
  fIsFluxWghtUpdatedRHC = false;
  //fIsTweakedFHCpRHC         = true;
  //fIsFluxWghtUpdatedFHCpRHC = false;
  //

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the banff covariance matrix
  //

  if(fFHCpRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeFHCpRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  // NEW FHC+RHC JOINT
  TVectorD & bestfit_values_fhc = *fBestFitValues;
  TVectorD & current_values_fhc = *fCurrentValues;
  TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  TVectorD & current_values_rhc = *fRHCCurrentValues;
  //

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
  {
     // Input vector is NULL.
     // This resets the BANFF params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakXSecOnly >> Reset FHCpRHC BANFF xsec params to best-fit values." << endl;
#endif
	// NEW FHC+RHC JOINT
        fIsTweaked    = false;
        fIsTweakedRHC = false;
        //fIsTweakedFHCpRHC = false;
	//
        
	current_values = bestfit_values;
	// NEW FHC+RHC JOINT
        current_values_fhc = bestfit_values_fhc;
        current_values_rhc = bestfit_values_rhc;
        //
        return;
     }
     else     
     {
#ifdef __show_tweak_mesg_debug__
       cout << "Exp/BANFF/TweakXSecOnly >> Forcing input FHCpRHC BANFF xsec parameters. " << endl;
#endif
       TVectorD & vec = *input_vec;
       int knumofxsecparamsinbanfffit = conv::kNumOfFHCpRHCParamsInBANFFfit - conv::kNumOfFHCpRHCFluxParamsInBANFFfit; 
       assert(vec.GetNoElements() == knumofxsecparamsinbanfffit);
       for (int i=0; i< knumofxsecparamsinbanfffit; i++){
	 current_values[i] = vec[i];
       }

       // NEW FHC+RHC JOINT
       this->SplitFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);
      //

       return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/TweakFluxOnly >> Generating FHCpRHC random twk dials (correlated)! " << endl;
#endif
     TMatrixD & L = *fFHCpRHCCovMtrxChoL;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);

     for (int i=0; i< conv::kNumOfFHCpRHCParamsInBANFFfit; i++){
       if(i>=conv::kNumOfFHCpRHCFluxParamsInBANFFfit) current_values[i] = vec[i]; // xsec params are after flux params
       else current_values[i] = bestfit_values[i];
     }

     // NEW FHC+RHC JOINT
     this->SplitFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);
     //
     
     return;
  }
}




//________________________________________________________________________________________
void BANFF::TweakSampleFluxOnly(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
  #if defined(__FULLY_CORRELATED__) 
  this->TweakFluxOnly(force_inputs,input_vec);
  this->SyncRHCToFHCSystematics();
  #else
  this->TweakFHCpRHCFluxOnly(force_inputs,input_vec);
  #endif
  //
#else
#ifdef __FHC_RUNNING__
  this->TweakFluxOnly(force_inputs,input_vec);
#endif
#ifdef __RHC_RUNNING__
  this->TweakRHCFluxOnly(force_inputs,input_vec);
#endif
#endif
}
//________________________________________________________________________________________
void BANFF::TweakFluxOnly(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/TweakFluxOnly >> Tweaking BANFF flux parameters" << endl;
#endif

  // Switch flags
  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fCovMtrxChoL == 0) 
    this->CholeskyDecomposeMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fBestFitValues;
  TVectorD & current_values = *fCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
  {
     // Input vector is NULL.
     // This resets the BANFF params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakFluxOnly >> Reset BANFF flux params to best-fit values." << endl;
#endif
        fIsTweaked = false;
        current_values = bestfit_values;
        return;
     }
     else     
     {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakFluxOnly >> Forcing input BANFF flux parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;
        assert(vec.GetNoElements() == conv::kNumOfFluxParamsInBANFFfit);
	for (int i=0; i< conv::kNumOfFluxParamsInBANFFfit; i++){
	  current_values[i] = vec[i];
	}
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/TweakFluxOnly >> Generating random twk dials (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
     for (int i=0; i< conv::kNumOfParamsInBANFFfit; i++){
       if(i<conv::kNumOfFluxParamsInBANFFfit) current_values[i] = vec[i];
       else current_values[i] = bestfit_values[i];
     }
     return;
  }

}
//________________________________________________________________________________________
void BANFF::TweakRHCFluxOnly(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/TweakFluxOnly >> Tweaking RHC BANFF flux parameters" << endl;
#endif

  // Switch flags
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdatedRHC = false;

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fRHCBestFitValues;
  TVectorD & current_values = *fRHCCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakFluxOnly >> Reset RHC BANFF flux params to best-fit values." << endl;
#endif
        fIsTweakedRHC = false;
        current_values = bestfit_values;
        return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/TweakFluxOnly >> Forcing RHC input BANFF flux parameters. " << endl;
#endif
	  TVectorD & vec = *input_vec;
	  assert(vec.GetNoElements() == conv::kNumOfRHCFluxParamsInBANFFfit);
	  for (int i=0; i< conv::kNumOfRHCFluxParamsInBANFFfit; i++){
	    current_values[i] = vec[i];
	  }
	  return;
	}
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/TweakFluxOnly >> Generating RHC random twk dials (correlated)! " << endl;
#endif
      TMatrixD & L = *fRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
      for (int i=0; i< conv::kNumOfRHCParamsInBANFFfit; i++){
	if(i<conv::kNumOfRHCFluxParamsInBANFFfit) current_values[i] = vec[i];
	else current_values[i] = bestfit_values[i];
      }
      return;
    }
}
//________________________________________________________________________________________





// NEW FLUX FHC+RHC

void BANFF::TweakFHCpRHCFluxOnly(
 bool       force_inputs, /* force the input params or generate random ones internally  */
 TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/TweakFluxOnly >> Tweaking FHC+RHC BANFF flux parameters" << endl;
#endif

  // Switch flags

  // NEW FHC+RHC JOINT
  fIsTweaked            = true;
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdated    = false;
  fIsFluxWghtUpdatedRHC = false;
  //fIsTweakedFHCpRHC         = true;
  //fIsFluxWghtUpdatedFHCpRHC = false;
  //

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //

  if(fFHCpRHCCovMtrxChoL == 0) 
    this->CholeskyDecomposeFHCpRHCMatrix();

  //
  // Get best-fit and current BANFF params
  //

  TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  // NEW FHC+RHC JOINT
  TVectorD & bestfit_values_fhc = *fBestFitValues;
  TVectorD & current_values_fhc = *fCurrentValues;
  TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  TVectorD & current_values_rhc = *fRHCCurrentValues;
  //

  //
  // Force input BANFF parameters
  //
  if(force_inputs) 
    {
      // Input vector is NULL.
      // This resets the BANFF params to their best-fit values.
      if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/TweakFluxOnly >> Reset FHC+RHC BANFF flux params to best-fit values." << endl;
#endif
        // NEW FHC+RHC JOINT
	fIsTweaked          = false;
        fIsTweakedRHC       = false;
        //fIsTweakedFHCpRHC = false;
        
	current_values = bestfit_values;
	// NEW FHC+RHC JOINT
        current_values_fhc = bestfit_values_fhc;
        current_values_rhc = bestfit_values_rhc;
        //
        return;
      }
      else     
	{
#ifdef __show_tweak_mesg_debug__
	  cout << "Exp/BANFF/TweakFluxOnly >> Forcing FHC+RHC input BANFF flux parameters. " << endl;
#endif
	  TVectorD & vec = *input_vec;
	  assert(vec.GetNoElements() == conv::kNumOfFHCpRHCFluxParamsInBANFFfit);
	  for (int i=0; i< conv::kNumOfFHCpRHCFluxParamsInBANFFfit; i++){
	    current_values[i] = vec[i];
	  }

	  // NEW FHC+RHC JOINT
	  this->SplitFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);
	  //

	  return;
	}
    }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
    {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/TweakFluxOnly >> Generating FHC+RHC random twk dials (correlated)! " << endl;
#endif
      TMatrixD & L = *fFHCpRHCCovMtrxChoL;
      TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values);
      for (int i=0; i< conv::kNumOfFHCpRHCParamsInBANFFfit; i++){
	if(i<conv::kNumOfFHCpRHCFluxParamsInBANFFfit) current_values[i] = vec[i];
	else current_values[i] = bestfit_values[i];
      }

      // NEW FHC+RHC JOINT
      this->SplitFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);
      //
      
      return;
    }
}
//

//________________________________________________________________________________________
void BANFF::CholeskyDecomposeMatrix()
{
#ifdef __show_cholesky_mesg_info__  
    cout << "Exp/BANFF/CholeskyDecomposeMatrix >> "
         << "Running Cholesky's decomposition method on the BANFF covariance matrix." << endl;
#endif

    TMatrixD cov(this->CovMatrix());

#ifdef __show_cholesky_mesg_debug__
    cout << "Exp/BANFF/CholeskyDecomposeMatrix >> BANFF covariance matrix: " << endl; 
    util::PrintCompact(cov,7,5);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_cholesky_mesg_debug__
    cout << "Exp/BANFF/CholeskyDecomposeMatrix >> Cholesky decomposition of BANFF covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,5);
#endif
}
//________________________________________________________________________________________
void BANFF::CholeskyDecomposeRHCMatrix()
{
#ifdef __show_cholesky_mesg_info__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> "
       << "Running Cholesky's decomposition method on the RHC BANFF covariance matrix." << endl;
#endif

  TMatrixD cov(this->RHCCovMatrix());
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(cov,7,5);
#endif

  TMatrixD L = util::CholeskyDecomposition(cov);
  fRHCCovMtrxChoL = new TMatrixD(L);
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> Cholesky decomposition of RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(*fRHCCovMtrxChoL,7,5);
#endif
}
//________________________________________________________________________________________



// NEW FLUX FHC+RHC
void BANFF::CholeskyDecomposeFHCpRHCMatrix()
{
#ifdef __show_cholesky_mesg_info__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> "
       << "Running Cholesky's decomposition method on the FHC+RHC BANFF covariance matrix." << endl;
#endif

  TMatrixD cov(this->FHCpRHCCovMatrix());
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> FHC+RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(cov,7,5);
#endif

  TMatrixD L = util::CholeskyDecomposition(cov);
  fFHCpRHCCovMtrxChoL = new TMatrixD(L);
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/CholeskyDecomposeMatrix >> Cholesky decomposition of FHC+RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(*fFHCpRHCCovMtrxChoL,7,5);
#endif
}
//

//________________________________________________________________________________________
void BANFF::RestoreNominalMC(void)
{
// A special tweak method

  cout << "Exp/BANFF >> ** Applying special tweak cancelling the effect of BANFF tuning "
       << "and restoring nominal MC..." << endl; 

  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //this is doing nothing at end, because NuFlux and NuXSec
  //are still reading bestfit values from fBestFitValues
  //we need to modify BestFitValue(p) 
  TVectorD & current_values = *fCurrentValues;
  //~ TVectorD & bestfit_values = *fBestFitValues;

  const int nfp = conv::kNumOfFluxParamsInBANFFfit;
  for(int i=0; i<nfp; i++) {
    current_values[i] = 1.;
  }

  const int nxsp = conv::kNumOfParamsInBANFFfit - kNumOfFluxParamsInBANFFfit;
  for (int i = 0 ; i < nxsp ; i++) {
    current_values[nfp+i] = (*fUntunedCentralValues)[i]; 
  }
}
//________________________________________________________________________________________
void BANFF::RestoreRHCNominalMC(void)
{
  // A special tweak method

  cout << "Exp/BANFF >> ** Applying special tweak cancelling the effect of RHC BANFF tuning "
       << "and restoring RHC nominal MC..." << endl; 

  fIsTweakedRHC         = true;
  fIsFluxWghtUpdatedRHC = false;

  //this is doing nothing at end, because NuFlux and NuXSec
  //are still reading bestfit values from fBestFitValues
  //we need to modify BestFitValue(p) 
  TVectorD & current_values = *fRHCCurrentValues;
  //~ TVectorD & bestfit_values = *fRHCBestFitValues;

  const int nfp = conv::kNumOfRHCFluxParamsInBANFFfit;
  for(int i=0; i<nfp; i++) {
    current_values[i] = 1.;
  }
  const int nxsp = conv::kNumOfRHCParamsInBANFFfit - kNumOfRHCFluxParamsInBANFFfit;
  for (int i = 0 ; i < nxsp ; i++) {
    current_values[nfp+i] = (*fRHCUntunedCentralValues)[i];
  }
}
//________________________________________________________________________________________

// NEW FHC+RHC JOINT
void BANFF::RestoreFHCpRHCNominalMC(void)
{
  // A special tweak method

  cout << "Exp/BANFF >> ** Applying special tweak cancelling the effect of FHC+RHC BANFF tuning "
       << "and restoring FHC+RHC nominal MC..." << endl; 

  // NEW FHC+RHC JOINT
  fIsTweaked            = true;
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdated    = false;
  fIsFluxWghtUpdatedRHC = false;
  //fIsTweakedFHCpRHC         = true;
  //fIsFluxWghtUpdatedFHCpRHC = false;
  //
  
  this->RestoreNominalMC();
  this->RestoreRHCNominalMC();
  
  //TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  //TVectorD & bestfit_values_fhc = *fBestFitValues;
  TVectorD & current_values_fhc = *fCurrentValues;
  //TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  TVectorD & current_values_rhc = *fRHCCurrentValues;
  this->MergeFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);    
}
//

//________________________________________________________________________________________
void BANFF::RestoreTunedMC(void)
{
// A special tweak method

  cout << "Exp/BANFF >> ** Applying special tweak applying the effect of BANFF tuning "
       << "and restoring tuned MC..." << endl; 

  fIsTweaked         = false;
  fIsFluxWghtUpdated = true;

  TVectorD & bestfit_values = *fBestFitValues;
  TVectorD & current_values = *fCurrentValues;
  current_values = bestfit_values;
}

//________________________________________________________________________________________
void BANFF::RestoreRHCTunedMC(void)
{
// A special tweak method

  cout << "Exp/BANFF >> ** Applying special tweak applying the effect of RHC BANFF tuning "
       << "and restoring RHC tuned MC..." << endl; 

  fIsTweakedRHC         = false;
  fIsFluxWghtUpdatedRHC = true;

  TVectorD & bestfit_values = *fRHCBestFitValues;
  TVectorD & current_values = *fRHCCurrentValues;
  current_values = bestfit_values;
}
//________________________________________________________________________________________

// NEW FHC+RHC JOINT

void BANFF::RestoreFHCpRHCTunedMC(void)
{
  // A special tweak method
  
  cout << "Exp/BANFF >> ** Applying special tweak cancelling the effect of FHC+RHC BANFF tuning "
       << "and restoring FHC+RHC nominal MC..." << endl; 
  
  // NEW FHC+RHC JOINT
  fIsTweaked            = true;
  fIsTweakedRHC         = true;
  fIsFluxWghtUpdated    = false;
  fIsFluxWghtUpdatedRHC = false;
  //fIsTweakedFHCpRHC         = true;
  //fIsFluxWghtUpdatedFHCpRHC = false;
  //

  this->RestoreTunedMC();
  this->RestoreRHCTunedMC();
  
  //TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  //TVectorD & bestfit_values_fhc = *fBestFitValues;
  TVectorD & current_values_fhc = *fCurrentValues;
  //TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  TVectorD & current_values_rhc = *fRHCCurrentValues;
  this->MergeFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);    
}
//


//________________________________________________________________________________________
double BANFF::ChisqPenalty(bool /*all_syst*/) const
{
#if !defined(__FHC_RUNNING__) && !defined(__MULTIRING_FHC_RUNNING__)
  return 0;
#endif

  //Vector of current values & best fit values
  const TVectorD & bestfit_values   = *fBestFitValues;
  const TVectorD & current_values   = *fCurrentValues;
#ifdef __show_chisqpenalty_mesg_debug2__
    cout << "Exp/BANFF/ChisqPenalty >> BANFF parameters: " << endl;
    util::PrintCompact(current_values, bestfit_values, 12,6,"BANFF fit params");
#endif

  //covariance matrix
  const TMatrixD cov    (this->CovMatrix());
  const TMatrixD covInv (this->CovMatrixInv());

  //Vector with deviations of current values from  best fit values
  const int np = conv::kNumOfParamsInBANFFfit;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  for(int i = 0; i < np; i++){
    dev[0][i]  = current_values[i]-bestfit_values[i];
    devT[i][0] = current_values[i]-bestfit_values[i];
  }
  
#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in BANFF parameters: " << endl;
  util::PrintCompact(TMatrixDRow(dev,0), 12,6,"BANFF fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug2__
    cout << "Exp/BANFF/ChisqPenalty >> BANFF covariance matrix: " << endl;
    util::PrintCompact(cov,7,5);
    cout << "Exp/BANFF/ChisqPenalty >> Inverse BANFF covariance matrix: " << endl;
    util::PrintCompact(covInv,7,5);
#endif

    //Check the inverse:
#ifdef __show_chisqpenalty_mesg_debug2__
    TMatrixD prod = cov*covInv;
    cout << "Exp/BANFF/ChisqPenalty >> check of product cov*(cov)^-1 : "<< endl;
    util::PrintCompact(prod,7,5);
#endif

    TMatrixD chisqpenalty =  dev * covInv * devT;
    double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

#ifdef __show_chisqpenalty_mesg_debug__
    cout << "Exp/BANFF/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;
#endif

    if(chisq_penalty) return chisq_penalty;
    //}//all_syst

  return 0;
}
//________________________________________________________________________________________ 
double BANFF::RHCChisqPenalty(bool /*all_syst*/) const
{
#ifndef __RHC_RUNNING__
  return 0;
#endif
  //Vector of current values & best fit values
  const TVectorD & bestfit_values   = *fRHCBestFitValues;
  const TVectorD & current_values   = *fRHCCurrentValues;
#ifdef __show_chisqpenalty_mesg_debug2__
  cout << "Exp/BANFF/ChisqPenalty >> RHC BANFF parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"RHC BANFF fit params");
#endif

  //covariance matrix 
  const TMatrixD cov    (this->RHCCovMatrix());
  const TMatrixD covInv (this->RHCCovMatrixInv());

  //Vector with deviations of current values from  best fit values
  const int np = conv::kNumOfRHCParamsInBANFFfit;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  for(int i = 0; i < np; i++){
    dev[0][i]  = current_values[i]-bestfit_values[i];
    devT[i][0] = current_values[i]-bestfit_values[i];
  }

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in RHC BANFF parameters: " << endl;
  util::PrintCompact(TMatrixDRow(dev,0), 12,6,"RHC BANFF fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug2__
  cout << "Exp/BANFF/ChisqPenalty >> RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF/ChisqPenalty >> RHC inverse BANFF covariance matrix: " << endl;
  util::PrintCompact(covInv,7,5);
#endif

  //Check the inverse:
#ifdef __show_chisqpenalty_mesg_debug2__
  TMatrixD prod = cov*covInv;
  cout << "Exp/BANFF/ChisqPenalty >> check of RHC product cov*(cov)^-1 : "<< endl;
  util::PrintCompact(prod,7,5);
#endif

  TMatrixD chisqpenalty =  dev * covInv * devT;
  double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/BANFF/ChisqPenalty >> RHC ChisqPenalty = " << chisq_penalty << endl;
#endif

  if(chisq_penalty) return chisq_penalty;
  //}//all_syst

  return 0;
}
//________________________________________________________________________________________ 





// NEW FLUX FHC+RHC
//double BANFF::FHCpRHCChisqPenalty(bool /*all_syst*/) const
double BANFF::FHCpRHCChisqPenalty(bool /*all_syst*/)
{
#if !(defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__))
  return 0;
#endif

  //Vector of current values & best fit values
  TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  
  // NEW FHC+RHC JOINT  
  const TVectorD & bestfit_values_fhc = *fBestFitValues;
  const TVectorD & current_values_fhc = *fCurrentValues;
  const TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  const TVectorD & current_values_rhc = *fRHCCurrentValues;
  // Merge FHC and RHC vectors into a single FHCpRHC vector 
  // Take into account FHC/RHC correlations
  this->MergeFHCpRHCvectors(bestfit_values,bestfit_values_fhc,bestfit_values_rhc);
  this->MergeFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);
  //

#ifdef __show_chisqpenalty_mesg_debug2__
  cout << "Exp/BANFF/ChisqPenalty >> FHC+RHC BANFF parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"FHC+RHC BANFF fit params");
#endif

  //covariance matrix 
  const TMatrixD cov    (this->FHCpRHCCovMatrix());
  const TMatrixD covInv (this->FHCpRHCCovMatrixInv());

  //Vector with deviations of current values from  best fit values
  const int np = conv::kNumOfFHCpRHCParamsInBANFFfit;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  for(int i = 0; i < np; i++){
    dev[0][i]  = current_values[i]-bestfit_values[i];
    devT[i][0] = current_values[i]-bestfit_values[i];
  }

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/ChisqPenalty >> deviations in FHC+RHC BANFF parameters: " << endl;
  util::PrintCompact(TMatrixDRow(dev,0), 12,6,"FHC+RHC BANFF fit params");
#endif

#ifdef __show_chisqpenalty_mesg_debug2__
  cout << "Exp/BANFF/ChisqPenalty >> FHC+RHC BANFF covariance matrix: " << endl;
  util::PrintCompact(cov,7,5);
  cout << "Exp/BANFF/ChisqPenalty >> FHC+RHC inverse BANFF covariance matrix: " << endl;
  util::PrintCompact(covInv,7,5);
#endif

  //Check the inverse:
#ifdef __show_chisqpenalty_mesg_debug2__
  TMatrixD prod = cov*covInv;
  cout << "Exp/BANFF/ChisqPenalty >> check of FHC+RHC product cov*(cov)^-1 : "<< endl;
  util::PrintCompact(prod,7,5);
#endif

  TMatrixD chisqpenalty =  dev * covInv * devT;
  double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

#ifdef __show_chisqpenalty_mesg_debug__
  cout << "Exp/BANFF/ChisqPenalty >> FHC+RHC ChisqPenalty = " << chisq_penalty << endl;
#endif

  if(chisq_penalty) return chisq_penalty;
  //}//all_syst

  return 0;
}

//________________________________________________________________________________________
void BANFF::PrintSystematicVariations(std::ostream & stream) const
{
  //assert(fCurrentValues->GetNoElements() == conv::kNumOfParamsInBANFFfit);
  const TVectorD & current_values = *fCurrentValues;
  const TVectorD & bestfit_values = *fBestFitValues;
  stream << "Exp/BANFF >> BANFF parameters: " << endl; 
  util::PrintCompact(current_values, bestfit_values, 12,6,"FHC BANFF fit params", stream);
}
//________________________________________________________________________________________
void BANFF::PrintRHCSystematicVariations(std::ostream & stream) const
{
  //assert(fCurrentValues->GetNoElements() == conv::kNumOfRHCParamsInBANFFfit);
  const TVectorD & current_values = *fRHCCurrentValues;
  const TVectorD & bestfit_values = *fRHCBestFitValues;
  stream << "Exp/BANFF >> RHC BANFF parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"RHC BANFF fit params", stream);
}
//________________________________________________________________________________________



// NEW FLUX FHC+RHC
void BANFF::PrintFHCpRHCSystematicVariations(std::ostream & stream) const
{
  //assert(fCurrentValues->GetNoElements() == conv::kNumOfRHCParamsInBANFFfit);
  const TVectorD & current_values = *fFHCpRHCCurrentValues;
  const TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  stream << "Exp/BANFF >> FHC+RHC BANFF parameters: " << endl;
  util::PrintCompact(current_values, bestfit_values, 12,6,"FHC+RHC BANFF fit params", stream);
}
//





// NEW FHC+RHC JOINT

//________________________________________________________________________________________
int BANFF::Param2Idx_FHCpRHC(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  return fParam2Idx_FHCpRHC[p];  
}
//________________________________________________________________________________________
int BANFF::Param2Idx_FHC(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  return fParam2Idx_FHC[p];  
}
int BANFF::Param2Idx_RHC(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  return fParam2Idx_RHC[p];  
}
//

//________________________________________________________________________________________
const TVectorD & BANFF::SampleBestFitValues (void) const
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) 
  return *fBestFitValues;
  #else
  return *fFHCpRHCBestFitValues;
  #endif
  //
#elif defined(__FHC_RUNNING__)
  return *fBestFitValues;
#elif defined(__RHC_RUNNING__)
  return *fRHCBestFitValues;
#else
  std::cerr<<"BANFF/SampleBestFitValues >> Cannot find right matrix."<<endl;
  exit(-1);
#endif
}
//________________________________________________________________________________________
const TMatrixD & BANFF::SampleCovMatrix  (void) const
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) 
  return *fCovMtrx;
  #else
  return *fFHCpRHCCovMtrx;
  #endif
  //
#elif defined(__FHC_RUNNING__)
  return *fCovMtrx;
#elif defined(__RHC_RUNNING__)
  return *fRHCCovMtrx;
#else
  std::cerr<<"BANFF/SampleCovMtrx >> Cannot find right matrix."<<endl;
  exit(-1);
#endif
}
//________________________________________________________________________________________
const TMatrixD & BANFF::SampleCovMatrixInv (void) const
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) 
  return *fCovMtrxInv;
  #else
  return *fFHCpRHCCovMtrxInv;
  #endif
  //
#elif defined(__FHC_RUNNING__)
  return *fCovMtrxInv;
#elif defined(__RHC_RUNNING__)
  return *fRHCCovMtrxInv;
#else
  std::cerr<<"BANFF/SampleCovMtrxInv >> Cannot find right matrix."<<endl;
  exit(-1);
#endif
}
//________________________________________________________________________________________




// NEW FHC+RHC JOINT

void BANFF::SplitFHCpRHCvectors(TVectorD bestfit_values_fhcprhc,TVectorD &bestfit_values_fhc,TVectorD &bestfit_values_rhc){
  
  // Fill FHC and RHC parameter arrays

  const int np = conv::kNumOfFHCpRHCParamsInBANFFfit;

  // Loop over all fhc+rhc parameters and fill fhc- and rhc-only vectors
  // The order of the parameters in the input vectors is: flux fhc, flux rhc, xsec fhc/rhc common
  int range_flux_fhc[2] = {0, conv::kNumOfFluxParamsInBANFFfit}; // flux fhc
  int range_flux_rhc[2] = {conv::kNumOfFluxParamsInBANFFfit,    conv::kNumOfFHCpRHCFluxParamsInBANFFfit}; // flux rhc 
  int range_xsec_fhc[2] = {conv::kNumOfFHCpRHCFluxParamsInBANFFfit, conv::kNumOfFHCpRHCParamsInBANFFfit}; // xsec fhc
  int range_xsec_rhc[2] = {conv::kNumOfFHCpRHCFluxParamsInBANFFfit, conv::kNumOfFHCpRHCParamsInBANFFfit}; // xsec rhc
  
#ifdef __show_init_mesg_debug__
  cout << "range flux fhc = " <<  range_flux_fhc[0] << " - " <<  range_flux_fhc[1] << endl;
  cout << "range flux rhc = " <<  range_flux_rhc[0] << " - " <<  range_flux_rhc[1] << endl;
  cout << "range xsec fhc = " <<  range_xsec_fhc[0] << " - " <<  range_xsec_fhc[1] << endl;
  cout << "range xsec rhc = " <<  range_xsec_rhc[0] << " - " <<  range_xsec_rhc[1] << endl;
#endif

  int count_fhc=0;
  int count_rhc=0;
  for(int ipar=0;ipar<np;ipar++){
    
    // use only if statements because vectors can share same parameters (e.g. xsec)

    // fhc
    if( ipar >= range_flux_fhc[0] && ipar < range_flux_fhc[1] ){     
      bestfit_values_fhc[ count_fhc ] = bestfit_values_fhcprhc[ipar];
      count_fhc++;
    }
    if( ipar >= range_xsec_fhc[0] && ipar < range_xsec_fhc[1] ){
      bestfit_values_fhc[ count_fhc ] = bestfit_values_fhcprhc[ipar];
      count_fhc++;
    }
    // rhc
    if( ipar >= range_flux_rhc[0] && ipar < range_flux_rhc[1] ){
      bestfit_values_rhc[ count_rhc ] = bestfit_values_fhcprhc[ipar];
      count_rhc++;
    }    
    if( ipar >= range_xsec_rhc[0] && ipar < range_xsec_rhc[1] ){
      bestfit_values_rhc[ count_rhc ] = bestfit_values_fhcprhc[ipar];
      count_rhc++;
    }    
  } 
 
  //bestfit_values_fhc.Print();
  //bestfit_values_rhc.Print();
  //bestfit_values_fhcprhc.Print();

  return;
}


void BANFF::MergeFHCpRHCvectors(TVectorD &bestfit_values_fhcprhc,TVectorD bestfit_values_fhc,TVectorD bestfit_values_rhc)
{
  
  // Fill FHC and RHC parameter arrays
  
  //const int np = conv::kNumOfFHCpRHCParamsInBANFFfit;
  const int np_fhc = conv::kNumOfParamsInBANFFfit;
  const int np_rhc = conv::kNumOfRHCParamsInBANFFfit;
  
  // Loop over all fhc+rhc parameters and fill fhc- and rhc-only vectors
  // The order of the parameters in the input vectors is: flux fhc, flux rhc, xsec fhc/rhc common
  int range_flux_fhc[2] = {0, conv::kNumOfFluxParamsInBANFFfit}; // flux fhc
  int range_flux_rhc[2] = {0, conv::kNumOfRHCFluxParamsInBANFFfit}; // flux rhc 
  int range_xsec_fhc[2] = {conv::kNumOfFluxParamsInBANFFfit,    conv::kNumOfParamsInBANFFfit}; // xsec fhc

#ifdef __show_init_mesg_debug__
  int range_xsec_rhc[2] = {conv::kNumOfRHCFluxParamsInBANFFfit, conv::kNumOfRHCParamsInBANFFfit}; // xsec rhc
  
  cout << "range flux fhc = " <<  range_flux_fhc[0] << " - " <<  range_flux_fhc[1] << endl;
  cout << "range flux rhc = " <<  range_flux_rhc[0] << " - " <<  range_flux_rhc[1] << endl;
  cout << "range xsec fhc = " <<  range_xsec_fhc[0] << " - " <<  range_xsec_fhc[1] << endl;
  cout << "range xsec rhc = " <<  range_xsec_rhc[0] << " - " <<  range_xsec_rhc[1] << endl;
#endif 

  int count=0;
  
  // fhc flux
  for(int ipar=0;ipar<np_fhc;ipar++){
    if( ipar >= range_flux_fhc[0] && ipar < range_flux_fhc[1] ){     
      bestfit_values_fhcprhc[count] = bestfit_values_fhc[ ipar ]; 
      count++;
    }
  }

  // rhc flux
  for(int ipar=0;ipar<np_rhc;ipar++){
    if( ipar >= range_flux_rhc[0] && ipar < range_flux_rhc[1] ){     
      bestfit_values_fhcprhc[count] = bestfit_values_rhc[ ipar ]; 
      count++;
    }
  }
  
  // fhc xsec (same as rhc xsec)
  for(int ipar=0;ipar<np_fhc;ipar++){
    if( ipar >= range_xsec_fhc[0] && ipar < range_xsec_fhc[1] ){     
      bestfit_values_fhcprhc[count] = bestfit_values_fhc[ ipar ]; 
      count++;
    }
  }

  // bestfit_values_fhc.Print();
  // bestfit_values_rhc.Print();
  // bestfit_values_fhcprhc.Print();


  return;
}

void BANFF::UpdateFHCpRHCCurrentValues(){
  // Update best-fit and current values
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  TVectorD & current_values_fhc = *fCurrentValues;
  TVectorD & current_values_rhc = *fRHCCurrentValues;
  this->MergeFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);    
}

void BANFF::UpdateFHCpRHCBestFitValues(){
  // Update best-fit and current values
  TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & bestfit_values_fhc = *fBestFitValues;
  TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  this->MergeFHCpRHCvectors(bestfit_values,bestfit_values_fhc,bestfit_values_rhc);    
}

void BANFF::UpdateFHCandRHConlyCurrentValues(){
  // Update best-fit and current values
  TVectorD & current_values = *fFHCpRHCCurrentValues;
  TVectorD & current_values_fhc = *fCurrentValues;
  TVectorD & current_values_rhc = *fRHCCurrentValues;
  this->SplitFHCpRHCvectors(current_values,current_values_fhc,current_values_rhc);   
}

void BANFF::UpdateFHCandRHConlyBestFitValues(){
  // Update best-fit and current values
  TVectorD & bestfit_values = *fFHCpRHCBestFitValues;
  TVectorD & bestfit_values_fhc = *fBestFitValues;
  TVectorD & bestfit_values_rhc = *fRHCBestFitValues;
  this->SplitFHCpRHCvectors(bestfit_values,bestfit_values_fhc,bestfit_values_rhc);    
}


//
