//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance oscillation analysis
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  //
  PhysParamList fparm(2);
  fparm.AddParam(0, t2kana::yr1numu::conv::kPhypSinsq2Theta23, 0.9,   0.0,  1.0, 0.02, false, "sinsq2theta23");
  fparm.AddParam(1, t2kana::yr1numu::conv::kPhypDm23sq,        2E-3, 5E-4, 5E-2, 1E-5, false, "dmsq23");

  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  SuperkObservation skobs;
//  skobs.Read("inputs/skdata.root);

  OscillationFitter::FitInputs fit_inputs(skobs,skexp,fparm);

  OscillationFitter fitter(fit_inputs,"minuit");

  fitter.DoTheFit();
}
