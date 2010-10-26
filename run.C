{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libMinuit");   

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  PhysParamList fparm(2);
  fparm.AddParam(0, t2kana::yr1numu::conv::kPhypSinsq2Theta23, 0.9,   0.0,  1.0, 0.02, false, "sinsq2theta23");
  fparm.AddParam(1, t2kana::yr1numu::conv::kPhypDm23sq,        2E-3, 5E-4, 5E-2, 1E-5, false, "dmsq23");

  SuperkExpectation skexp;
//  skexp.GetND280Inputs()->Read("inputs/nd280.cfg");
//  skexp.GetSuperkInputs()->Read("inputs/superk.cfg");
//  skexp.GetBeamInputs()->Read("inputs/beam.cfg");
//  skexp.GetGeneratorInputs()->Read("inputs/genie.cfg");

  SuperkObservation skobs;
//  skobs.Read("inputs/skdata.root);

  OscillationFitter::FitInputs fit_inputs(skobs,skexp,fparm);

  OscillationFitter fitter(fit_inputs,"minuit");

  fitter.DoTheFit();
}
