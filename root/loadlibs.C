{
  TString fmwk_loadlibs = ".x $VALOR_FMWK/root/loadlibs.C";
  gSystem->ExpandPathName(fmwk_loadlibs);
  gROOT->ProcessLine(fmwk_loadlibs.Data());
  
  TString baselibs = gSystem->GetDynamicPath();
  TString t2k_libs = "$T2KOSCANA/lib:";
  t2k_libs += baselibs;
  
  gSystem->SetDynamicPath(t2k_libs.Data());
  
  gSystem->Load("libT2KNuOsc.so");
}
