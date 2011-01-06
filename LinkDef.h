#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace t2kana;
#pragma link C++ namespace t2kana::numu10a;
//#pragma link C++ namespace t2kana::numu10a::conv;
#pragma link C++ namespace t2kana::numu10a::util;

#pragma link C++ enum  t2kana::numu10a::conv::PhysParamId_t;
#pragma link C++ enum  t2kana::numu10a::conv::ConfLevel_t;
#pragma link C++ enum  t2kana::numu10a::conv::SuperkSample_t;
#pragma link C++ enum  t2kana::numu10a::conv::SuperkMode_t;

#pragma link C++ class t2kana::numu10a::SuperkObservation;
#pragma link C++ class t2kana::numu10a::SuperkExpectation;
#pragma link C++ class t2kana::numu10a::SuperkExpectation::ND280Inputs;
#pragma link C++ class t2kana::numu10a::SuperkExpectation::SuperkInputs;
#pragma link C++ class t2kana::numu10a::SuperkExpectation::BeamInputs;
#pragma link C++ class t2kana::numu10a::SuperkExpectation::GeneratorInputs;
#pragma link C++ class t2kana::numu10a::PhysParamList;
#pragma link C++ class t2kana::numu10a::OscillationFitter;
#pragma link C++ class t2kana::numu10a::OscillationFitter::FitInputs;
#pragma link C++ class t2kana::numu10a::DecayFitter;
#pragma link C++ class t2kana::numu10a::DecayFitter::FitInputs;
#pragma link C++ class t2kana::numu10a::DecohFitter;
#pragma link C++ class t2kana::numu10a::DecohFitter::FitInputs;

#endif
