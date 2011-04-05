{
  TString t2koscana_topdir = gSystem->Getenv("T2KOSCANA");

  TString mp = gROOT->GetMacroPath();
  TString ip;

  TString t2koscana_inc = t2koscana_topdir + "/fitter";

  const char* p = t2koscana_inc.Data();
  if (p) {
    mp += ":";
    mp += p;
    ip += " -I";
    ip += p;
  }

  mp += ":/usr/local/include/";
  ip += "  -I/usr/local/include/";

  gROOT->SetMacroPath(mp.Data());
  gSystem->SetIncludePath(ip);

  // additions to .include must be done individually or CINT will
  // try to quote all the spaces as a single path

  TString dip = ".include ";
  dip += t2koscana_inc.Data();
  gROOT->ProcessLine(dip.Data());

  dip = ".include /usr/local/include/";
  gROOT->ProcessLine(dip.Data());
}

