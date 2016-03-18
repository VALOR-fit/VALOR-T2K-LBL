{
  TString t2koscana_topdir = gSystem->Getenv("T2KOSCANA");
  TString valorfmwk_topdir = gSystem->Getenv("VALOR_FMWK");

  TString mp = gROOT->GetMacroPath();
  TString ip;

  TString t2koscana_inc = t2koscana_topdir + "/fitter";
  TString valorfmwk_inc = valorfmwk_topdir + "/src";

  const char* p = t2koscana_inc.Data();
  if (p) {
    mp += ":";
    mp += p;
    ip += " -I";
    ip += p;
  }
  const char* q = valorfmwk_inc.Data();
  if (q) {
    mp += ":";
    mp += q;
    ip += " -I";
    ip += q;
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

  dip = ".include ";
  dip += valorfmwk_inc.Data();
  gROOT->ProcessLine(dip.Data());

  dip = ".include /usr/local/include/";
  gROOT->ProcessLine(dip.Data());
}

