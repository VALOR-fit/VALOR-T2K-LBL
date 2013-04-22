#ifndef _t2kana_fitter_h_
#define _t2kana_fitter_h_

#include "FitterABC.h"
#include "FitterGSL.h"
#include "FitterMinuit.h"
#include "FitterMInuit2.h"
#include "FitInputs.h"

namespace t2kana {

//
// Option:
//   - minuit
//   - minuit2
//   - SimAnneal
//   - LevMarq
//   - ConjugateFR
//   - ConjugatePR
//   - VectorBFGS
//   - SteepestDescent
//

FitterABC * GetFitter (FitInputs * inputs, Option_t * option);

} // t2kana

#endif

