#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace t2kana;
#pragma link C++ namespace t2kana::conv;
#pragma link C++ namespace t2kana::util;

#pragma link C++ enum  t2kana::conv::PhysParamId_t;
#pragma link C++ enum  t2kana::conv::ConfLevel_t;
#pragma link C++ enum  t2kana::conv::SuperkSample_t;
#pragma link C++ enum  t2kana::conv::SuperkMode_t;

#pragma link C++ class t2kana::SuperkObservation;
#pragma link C++ class t2kana::SuperkExpectation;
#pragma link C++ class t2kana::SuperkMC;
#pragma link C++ class t2kana::ND280Measurement;
#pragma link C++ class t2kana::NuFlux;
#pragma link C++ class t2kana::NuXSec;
#pragma link C++ class t2kana::NuMuDisappFitter;
#pragma link C++ class t2kana::NuMuDecayFitter;
#pragma link C++ class t2kana::NuMuDecohFitter;
#pragma link C++ class t2kana::FitterABC;
#pragma link C++ class t2kana::FitParamList;
#pragma link C++ class t2kana::FitInputs;
#pragma link C++ class t2kana::PhysHypothesis;

#endif
