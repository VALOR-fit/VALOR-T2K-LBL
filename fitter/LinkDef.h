#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace valor;
#pragma link C++ namespace valor::conv;
#pragma link C++ namespace valor::util;

#pragma link C++ enum  valor::conv::PhysParamId_t;
#pragma link C++ enum  valor::conv::ConfLevel_t;
#pragma link C++ enum  valor::conv::SuperkSample_t;
#pragma link C++ enum  valor::conv::SuperkMode_t;


#pragma link C++ class valor::util::ArrIdx;
#pragma link C++ class valor::util::UtlMakeFakeData;

#pragma link C++ class valor::T2KPhysicalParameterLimits;
#pragma link C++ class valor::T2KOpt;

#pragma link C++ class valor::SuperkObservation;
#pragma link C++ class valor::SuperkExpectation;
#pragma link C++ class valor::SuperkDetectorFSI;
#pragma link C++ class valor::SuperkMC+;
#pragma link C++ class valor::BANFF;
#pragma link C++ class valor::NuXSec+;
#pragma link C++ class valor::NuFlux+;
#pragma link C++ class valor::FitMarginalization+;

#endif
