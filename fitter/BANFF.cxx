//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis) 
 see code in the corresponding SVN branch.
 
 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__  
//#define __show_init_mesg_debug__   
//#define __show_init_mesg_eigen_debug__
//#define __show_init_mesg_eigen_debug2__
#define __show_cholesky_mesg_info__
//#define __show_cholesky_mesg_debug__
//#define __show_tweak_mesg_debug__
//#define __show_update_correlated_mesg_debug__
//#define __show_chisqpenalty_mesg_debug__
//#define __show_chisqpenalty_mesg_debug2__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <TRandom3.h>
#include <TMath.h>

#include "BANFF.h" 
#include "Conventions.h" 
#include "Utilities.h" 

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

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
  fBestFitValues      = 0; 
  fCurrentValues      = 0; 
  fUncorrelatedCurrentValues = 0;
  fEigenValues        = 0;
  fEigenVectors       = 0;
  fEigenVectorsT      = 0;
  fEigenCovMtrx       = 0;
  fEigenBestFitValues = 0;
  fEigenCurrentValues = 0;
  fCovMtrx            = 0; 
  fCovMtrxChoL        = 0; 
  fIsTweaked          = false;
  fIsFluxWghtUpdated  = true;
}
//________________________________________________________________________________________
void BANFF::CleanUp(void)
{
  if ( fBestFitValues             ) delete fBestFitValues;
  if ( fCurrentValues             ) delete fCurrentValues;
  if ( fUncorrelatedCurrentValues ) delete fUncorrelatedCurrentValues;
  if ( fCovMtrx                   ) delete fCovMtrx;
  if ( fCovMtrxChoL               ) delete fCovMtrxChoL;
  if ( fEigenValues               ) delete fEigenValues;
  if ( fEigenVectors              ) delete fEigenVectors;
  if ( fEigenVectorsT             ) delete fEigenVectorsT;
  if ( fEigenCovMtrx              ) delete fEigenCovMtrx;
  if ( fEigenBestFitValues        ) delete fEigenBestFitValues;
  if ( fEigenCurrentValues        ) delete fEigenCurrentValues;
}
//________________________________________________________________________________________
bool BANFF::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read BANFF inputs to oscillation analysis
//
#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Reading BANFF inputs from XML configuration file." << endl;
#endif

  const int np = conv::kNumOfParamsInBANFFfit;

  // Get BANFF fit central values
  fBestFitValues = util::XML2TVectorD(xml,root_node,"banff/best_fit_values");
  if(!fBestFitValues) {
    cerr << "Exp/BANFF >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }

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
  fUncorrelatedCurrentValues = new TVectorD(bestfit);
  TVectorD & current_uncorrelated_values = *fUncorrelatedCurrentValues;
  for(int ip = 0; ip < np; ip++) {
    current_uncorrelated_values[ip] = 0;
  }

  // Get BANFF fit covariance matrix
  fCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/covariance_matrix");
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
  util::PrintCompact(cov,7,4);
  cout << "Exp/BANFF >> Inverse BANFF covariance matrix: " << endl;
  util::PrintCompact(*fCovMtrxInv,7,4);
#endif

  //EigenDecompose the BANFF fit covariance matrix
  TVectorD eigen_values(np);
  TMatrixD eigen_vectors = cov.EigenVectors(eigen_values);
  TMatrixD eigen_vectorsT(eigen_vectors);
  eigen_vectorsT.T();
  TMatrixDSym eigen_covariance(np);
  for(int i = 0; i < np; i++) {
    for(int j = 0; j < np; j++) {
      if(i==j)
	eigen_covariance(i,j) = eigen_values(i);
      else
	eigen_covariance(i,j) = 0.;
    }
  }
  TVectorD eigen_bestfit_values(bestfit);
  eigen_bestfit_values *= eigen_vectorsT;

  fEigenValues           = new TVectorD(eigen_values);
  fEigenVectors          = new TMatrixD(eigen_vectors);
  fEigenVectorsT         = new TMatrixD(eigen_vectorsT);
  fEigenCovMtrx          = new TMatrixDSym(eigen_covariance);
  fEigenBestFitValues    = new TVectorD(eigen_bestfit_values);
  //initialise current to bestfit
  fEigenCurrentValues    = new TVectorD(eigen_bestfit_values);

#ifdef __show_init_mesg_eigen_debug__
  eigen_values.Print();
  eigen_vectors.Print();
  eigen_vectorsT.Print();
  eigen_covariance.Print();
  eigen_bestfit_values.Print();
#endif

#ifdef __show_init_mesg_eigen_debug2__
  TVectorD & current_values = *fCurrentValues;
  current_values.Print();
  TVectorD & current_eigen_values = *fEigenCurrentValues;
  TVectorD vec(current_eigen_values);
  vec = eigen_vectors * current_eigen_values;
  vec.Print();
  exit(1);
#endif

#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Done reading BANFF inputs!" << endl;
#endif

  return true;
}
//________________________________________________________________________________________
double BANFF::BestFitValue(PhysParamId_t p) const
{
  const TVectorD & bestfit = this->BestFitValues();
  int idx = this->Param2Idx(p);
  return bestfit[idx];
}
//________________________________________________________________________________________
double BANFF::CurrentValue(PhysParamId_t p) const
{
  const TVectorD & current_values = *fCurrentValues;
  int idx = this->Param2Idx(p);
  return current_values[idx];
}
//________________________________________________________________________________________
double BANFF::OneSigmaError(PhysParamId_t p) const
{
  const TMatrixD & cov = *fCovMtrx;
  int idx = this->Param2Idx(p);
  return TMath::Sqrt(cov[idx][idx]);
}
//________________________________________________________________________________________
double BANFF::CurrentUncorrelatedValue(PhysParamId_t p) const
{
  const TVectorD & current_uncorrelated_values = *fUncorrelatedCurrentValues;
  int idx = this->UnParam2Idx(p);
  return current_uncorrelated_values[idx];
}
//________________________________________________________________________________________
double BANFF::CurrentEigenDecompValue(PhysParamId_t p) const
{
  const TVectorD & current_eigendecomp_values = *fEigenCurrentValues;
  int idx = this->EDParam2Idx(p);
  return current_eigendecomp_values[idx];
}
//________________________________________________________________________________________
double BANFF::BestFitEigenDecompValue(PhysParamId_t p) const
{
  const TVectorD & bestfit_eigendecomp_values = *fEigenBestFitValues;
  int idx = this->EDParam2Idx(p);
  return bestfit_eigendecomp_values[idx];
}
//________________________________________________________________________________________
double BANFF::EigenDecompEigenValue(PhysParamId_t p) const
{
  const TVectorD & eigen_values = *fEigenValues;
  int idx = this->EDParam2Idx(p);
  return eigen_values[idx];
}
//________________________________________________________________________________________
void BANFF::Tweak(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec,    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
  int        mode          /* tweak the 1: correlated (default) or 2: uncorrelated vector or 3: eigendecomp vector */
)
{

  /*
    WARNING. Do not mix working in correlated & uncorrelated parameter
    bases! When the correlated parameters are updated (BANFF::Tweak()),
    the uncorrelated parameters are not updated. If
    BANFF::UpdateCorrelated() (only available from BANFF::Tweak()) is then
    called, the result of the correlated tweak is lost.
    The same is true for eigendecomp!
  */

  if(mode == 1) {
    this->TweakCorrelated(force_inputs, input_vec);
  }
  else if(mode == 2) {
    this->TweakUncorrelated(force_inputs, input_vec);
    this->UpdateCorrelatedCholesky();
  }
  else if(mode == 3) {
    this->TweakEigenDecomp(force_inputs, input_vec);
    this->UpdateCorrelatedEigenDecomp();
  }
}
//________________________________________________________________________________________
void BANFF::TweakUncorrelated(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec    /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)
{
#ifdef __show_tweak_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking uncorrelated BANFF fit parameters" << endl;
#endif

  // Switch flags
  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //
  // Get current BANFF params
  //

  TVectorD & current_correlated_values = *fUncorrelatedCurrentValues;

  //
  // Force input BANFF parameters
  //
  if(force_inputs) {
    // Input vector is NULL.
    // This resets the BANFF params to their best-fit values.
    if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/Tweak >> Reset BANFF params to best-fit values." << endl;
#endif
      fIsTweaked = false;
      for(int ip = 0; ip < current_correlated_values.GetNrows(); ip++) {
	current_correlated_values[ip] = 0;
      }
      return;
    }
    else {
#ifdef __show_tweak_mesg_debug__
      cout << "Exp/BANFF/Tweak >> Forcing input uncorrelated BANFF parameters. " << endl;
#endif
      TVectorD & vec = *input_vec;
      assert(vec.GetNoElements() == conv::kNumOfParamsInBANFFfit);
      current_correlated_values = vec;
      return;
    }
  }//force_inputs?

  else{ 
#ifdef __show_tweak_mesg_debug__
    cout << "Exp/BANFF/Tweak >> Generating random twk dials (uncorrelated)!" << endl;
#endif

    TRandom3 rand(0);
    for(int ip = 0; ip < current_correlated_values.GetNrows(); ip++) {
      current_correlated_values[ip] = rand.Gaus(0,1);
    }
    return;
  }

}
//________________________________________________________________________________________
void BANFF::TweakCorrelated(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec     /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)    
{ 
#ifdef __show_cholesky_mesg_debug__
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
void BANFF::TweakEigenDecomp(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec     /* input vector of BANFF params - Set to 0 to reset to best-fit ones */
)    
{ 
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking eigendecomposition BANFF fit parameters" << endl;
#endif

  // Switch flags
  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  //
  // Get best-fit and current params
  //

  TVectorD & bestfit_eigen_values = *fEigenBestFitValues;
  TVectorD & current_eigen_values = *fEigenCurrentValues;

  //
  // Force input
  //
  if(force_inputs) 
  {
     // Input vector is NULL.
     // This resets the BANFF params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Reset eigendecomposition BANFF params to best-fit values." << endl;
#endif
        fIsTweaked = false;
        current_eigen_values = bestfit_eigen_values;
        return;
     }
     else     
     {
#ifdef __show_tweak_mesg_debug__
        cout << "Exp/BANFF/Tweak >> Forcing input eigendecomposition BANFF parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;
        assert(vec.GetNoElements() == conv::kNumOfParamsInBANFFfit);
        current_eigen_values = vec;
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate eigendecomposition tweaking dials
  //
  else
  {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/Tweak >> Generating random twk dials (eigendecomposition)! " << endl;
#endif
     TVectorD & eigen_values  = *fEigenValues;
     const int np = conv::kNumOfParamsInBANFFfit;
     TVectorD vec(np);
     TRandom3 rand(0);
     for(int ip = 0; ip < np; ip++) {
       vec(ip) = rand.Gaus(bestfit_eigen_values(ip), TMath::Sqrt(eigen_values(ip)));
     }
     current_eigen_values = vec;
     return;
  }

}
//________________________________________________________________________________________
void BANFF::UpdateCorrelatedCholesky()
{
#ifdef __show_update_correlated_mesg_debug__
  cout << "Exp/BANFF/UpdateCorrelatedCholesky >> Updating BANFF current values, using current uncorrelated parameters" << endl;
#endif

  //
  // Get best-fit & current BANFF params & decomposed cov matrix
  //
  TVectorD & bestfit_values = *fBestFitValues;
  TVectorD & current_values = *fCurrentValues;
  TVectorD & current_correlated_values = *fUncorrelatedCurrentValues;
  if(fCovMtrxChoL == 0)
    this->CholeskyDecomposeMatrix();
  TMatrixD & L = *fCovMtrxChoL;

  // Update the correlated current values
  TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,bestfit_values,current_correlated_values);
  current_values = vec;
  return;
}
//________________________________________________________________________________________
void BANFF::UpdateCorrelatedEigenDecomp()
{
#ifdef __show_update_correlated_mesg_debug__
  cout << "Exp/BANFF/UpdateCorrelatedEigenDecomp >> Updating BANFF current values, using current eigendecomposition parameters" << endl;
#endif

  //
  // Get current parameters, and transposed eigenvector matrix
  //
  TVectorD & current_values       = *fCurrentValues;
  TVectorD & current_eigen_values = *fEigenCurrentValues;
  TMatrixD & eigen_vectors        = *fEigenVectors;

  // Update the correlated current values
  TVectorD vec(current_eigen_values);
  vec = eigen_vectors * current_eigen_values;
  current_values = vec;
  return;
}
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
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_cholesky_mesg_debug__
    cout << "Exp/BANFF/CholeskyDecomposeMatrix >> Cholesky decomposition of BANFF covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,4);
#endif
}
//________________________________________________________________________________________
void BANFF::RestoreNominalMC(void)
{
// A special tweak method

// Does not work when working in eigendecomposition mode!

  cout << "Exp/BANFF >> ** Applying special tweak cancelling the effect of BANFF tuning "
       << "and restoring nominal MC..." << endl; 

  fIsTweaked         = true;
  fIsFluxWghtUpdated = false;

  TVectorD & current_values = *fCurrentValues;
  //TVectorD & current_uncorrelated_values = *fUncorrelatedCurrentValues;

  const int np = conv::kNumOfParamsInBANFFfit;
  for(int i=0; i<np; i++) {
    current_values[i] = 1.;
    //current_uncorrelated_values[i] = 0.;
  }
#if BANFF_BINS == 5
  current_values[np-1] = 0.; //mec!
#endif
}
//________________________________________________________________________________________
double BANFF::ChisqPenalty(bool all_syst, int mode) const
{
  double chisq_penalty = 0;
  if(mode == 1)
    chisq_penalty = ChisqPenaltyCorrelated(all_syst);
  else if(mode == 2)
    chisq_penalty = ChisqPenaltyUncorrelated(all_syst);
  else if(mode == 3)
    chisq_penalty = ChisqPenaltyEigenDecomp(all_syst);

  return chisq_penalty;
}
//________________________________________________________________________________________
double BANFF::ChisqPenaltyUncorrelated(bool all_syst) const
{
  double chisq_penalty = 0;
  if(!all_syst) {
    return 0;
    // Currently not fitting any uncorrelated BANFF parameters. If we do start, will use something like:
    //chisq_penalty += CurrentUncorrelatedValue(kNuisUnBANFF00) * CurrentUncorrelatedValue(kNuisUnBANFF00);
  }
  else {
    const TVectorD & current_correlated_values   = *fUncorrelatedCurrentValues;
    const int np = conv::kNumOfParamsInBANFFfit;
    for(int ip = 0; ip < np; ip++) {
      chisq_penalty += current_correlated_values[ip] * current_correlated_values[ip];
    }
  }

#ifdef __show_chisqpenalty_mesg_debug__
    cout << "Exp/BANFF/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;
#endif
    return chisq_penalty;
}
//________________________________________________________________________________________
double BANFF::ChisqPenaltyEigenDecomp(bool all_syst) const
{
  double chisq_penalty = 0;
  if(!all_syst) {
    return 0;
    // Currently not fitting any eigendecomposition BANFF parameters. If we do start, will use something like:
    //PhysParamId_t p = kNuisEDBANFF00;
    //double tmp = CurrentEigenDecompValue(p) - BestFitEigenDecompValue(p);
    //tmp *= tmp;
    //tmp /= EigenDecompEigenValue(p);
    //chisq_penalty += tmp;
  }
  else {
    const int np = conv::kNumOfParamsInBANFFfit;
    for(int ip = 0; ip < np; ip++) {
      PhysParamId_t p = PhysParamId_t ((int) kNuisEDBANFF00 + ip);
      double tmp = CurrentEigenDecompValue(p) - BestFitEigenDecompValue(p);
      tmp *= tmp;
      tmp /= EigenDecompEigenValue(p);
      chisq_penalty += tmp;
    }
  }

#ifdef __show_chisqpenalty_mesg_debug__
    cout << "Exp/BANFF/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;
#endif
    return chisq_penalty;
}
//________________________________________________________________________________________
double BANFF::ChisqPenaltyCorrelated(bool all_syst) const
{
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

  if(!all_syst) {
    //Vector with deviations of current values from  best fit values
    const int np = 2;
    TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
    TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
    int imaqe    = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFMAQE);
    int iccqeE0  = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFCCQEE0);
    dev [0][0] = current_values[imaqe]   - bestfit_values[imaqe];
    devT[0][0] = current_values[imaqe]   - bestfit_values[imaqe];
    dev [0][1] = current_values[iccqeE0] - bestfit_values[iccqeE0];
    devT[1][0] = current_values[iccqeE0] - bestfit_values[iccqeE0];

#ifdef __show_chisqpenalty_mesg_debug2__
    cout << "Exp/ChisqPenalty >> deviations in BANFF parameters: " << endl;
    util::PrintCompact(TMatrixDRow(dev,0), 12,6,"BANFF selected fit params");
#endif
      
    //fill sub cov matrix
    TMatrixD covSub(np, np);
    covSub[0][0] = cov[imaqe]   [imaqe];
    covSub[1][1] = cov[iccqeE0] [iccqeE0];
    covSub[0][1] = cov[imaqe]   [iccqeE0];
    covSub[1][0] = cov[iccqeE0] [imaqe];
    //invert the sub matrix
    double det;
    TMatrixD covInvSub(covSub);
    covInvSub.Invert(&det);

#ifdef __show_chisqpenalty_mesg_debug__
    cout << "Exp/BANFF/ChisqPenalty >> BANFF covariance sub matrix : " << endl;
    util::PrintCompact(covSub,7,4);
    cout << "Exp/BANFF/ChisqPenalty >> Inverse BANFF covariance sub matrix : " << endl;
    util::PrintCompact(covInvSub,7,4);
#endif

#ifdef __show_chisqpenalty_mesg_debug2__
    //check the inverse
    TMatrixD prod = covSub * covInvSub;
    cout << "Exp/BANFF/ChisqPenalty >> check of product covsub*(covsub)^-1 : "<< endl;
    util::PrintCompact(prod,7,4);
#endif

    TMatrixD chisqpenalty =  dev * covInvSub * devT;
    double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number
    return chisq_penalty;
  }

  else {
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
    util::PrintCompact(cov,7,4);
    cout << "Exp/BANFF/ChisqPenalty >> Inverse BANFF covariance matrix: " << endl;
    util::PrintCompact(covInv,7,4);
#endif

    //Check the inverse:
#ifdef __show_chisqpenalty_mesg_debug2__
    TMatrixD prod = cov*covInv;
    cout << "Exp/BANFF/ChisqPenalty >> check of product cov*(cov)^-1 : "<< endl;
    util::PrintCompact(prod,7,4);
#endif

    TMatrixD chisqpenalty =  dev * covInv * devT;
    double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number

#ifdef __show_chisqpenalty_mesg_debug__
    cout << "Exp/BANFF/ChisqPenalty >> ChisqPenalty = " << chisq_penalty << endl;
#endif

    if(chisq_penalty) return chisq_penalty;
  }

  return 0;
}
//________________________________________________________________________________________
void BANFF::PrintSystematicVariations(ostream & stream) const
{
  //assert(fCurrentValues->GetNoElements() == conv::kNumOfParamsInBANFFfit);

  const TVectorD & current_values = *fCurrentValues;
  const TVectorD & bestfit_values = *fBestFitValues;

  stream << "Exp/BANFF >> BANFF parameters: " << endl; 
  util::PrintCompact(current_values, bestfit_values, 12,6,"BANFF fit params", stream);
}
//________________________________________________________________________________________
int BANFF::Param2Idx(PhysParamId_t p) const
{
  int iparam = (int) p;
  int first  = (int) kNuisBANFF00;
#if BANFF_BINS <= 2
  int last   = (int) kNuisBANFFNC1pi0;
#elif BANFF_BINS == 5
  int last   = (int) kNuisBANFFMEC;
#elif BANFF_BINS >= 3
  int last   = (int) kNuisBANFFCC1piE1;
#endif
  if(! (iparam >= first && iparam <= last) )
    cout << AsString(p) << endl;
  assert(iparam >= first && iparam <= last);
  int idx = iparam-first;
  return idx;  
}
//________________________________________________________________________________________
int BANFF::UnParam2Idx(PhysParamId_t p) const
{
  int iparam = (int) p;
  int first  = (int) kNuisUnBANFF00;
  int last   = (int) kNuisUnBANFF00 + kNumOfParamsInBANFFfit;
  if(! (iparam >= first && iparam <= last) )
    cout << AsString(p) << endl;
  assert(iparam >= first && iparam <= last);
  int idx = iparam-first;
  return idx;  
}
//________________________________________________________________________________________
int BANFF::EDParam2Idx(PhysParamId_t p) const
{
  int iparam = (int) p;
  int first  = (int) kNuisEDBANFF00;
  int last   = (int) kNuisEDBANFF00 + kNumOfParamsInBANFFfit;
  if(! (iparam >= first && iparam <= last) )
    cout << AsString(p) << endl;
  assert(iparam >= first && iparam <= last);
  int idx = iparam-first;
  return idx;  
}
//________________________________________________________________________________________
