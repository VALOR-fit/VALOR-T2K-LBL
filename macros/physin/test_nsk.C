//
// ------------------------------------------------------------------------------------------------
// Print out the number of events predicted at SuperK (NSK)
// Also print out the number of events for certain categories. Using the same categories as
// table 16 in TN52 - which gives the results for run1+2 nue analysis B
// Update: now also produces (latex) table with number of events w & w/o BANFF tuning
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

//#include "Conventions.h"
#include <iomanip.h>

//structure to store results
struct nsk_t {
  double nsk = 0, nskmc = 0, nskmucc = 0, nskmubarcc = 0, nskecc = 0, nskebarcc = 0, nsksigecc = 0, nsknc = 0;
};

using namespace std;
using namespace t2kana;

nsk_t nsk_integration(nsk_t n, SuperkExpectation *skexp, conv::SuperkSample_t s);
void nsk_print(nsk_t n, double input_dm23sq);

void test_nsk
(
 double input_mixing13    = -99,
 double input_mixing23    = -99,
 double input_mixing12    = -99,
 double input_dm23sq      = -99,
 double input_dm12sq      = -99,
 double input_deltacp     = -99,
 const string sample_str  = "1Rmu",
 const char * config_file = "analysis_configuration.xml"
)
{
  conv::SuperkSample_t s = util::SKSampleFromString(sample_str);
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // fit params
  PhysHypothesis ph(t2kana::conv::kHypOsc3F);

  // generate predictions
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;

  ph.SetCurrValue(t2kana::conv::kPhypMixing23, input_mixing23);
  ph.SetCurrValue(t2kana::conv::kPhypMixing13, input_mixing13);
  ph.SetCurrValue(t2kana::conv::kPhypMixing12, input_mixing12);
  ph.SetCurrValue(t2kana::conv::kPhypDm23sq,   input_dm23sq);
  ph.SetCurrValue(t2kana::conv::kPhypDm12sq,   input_dm12sq);
  ph.SetCurrValue(t2kana::conv::kPhypDeltaCP,  input_deltacp);

  nsk_t osc_banff, noosc_banff, osc_nobanff, noosc_nobanff;

  // generate predictions
  skexp.Generate(ph);
  osc_banff = nsk_integration(osc_banff, &skexp, s);

  //undo banff
  skexp.GetBANFF()->RestoreNominalMC();
  skexp.Generate(ph);
  osc_nobanff = nsk_integration(osc_nobanff, &skexp, s);

  //restore banff, undo osc
  TVectorD twkv(skexp.GetBANFF()->BestFitValues());
  skexp.GetBANFF()->Tweak(true,&twkv);
  ph.SetCurrValue(t2kana::conv::kPhypMixing23, 0);
  ph.SetCurrValue(t2kana::conv::kPhypMixing13, 0);
  ph.SetCurrValue(t2kana::conv::kPhypMixing12, 0);
  ph.SetCurrValue(t2kana::conv::kPhypDm23sq,   0);
  ph.SetCurrValue(t2kana::conv::kPhypDm12sq,   0);
  ph.SetCurrValue(t2kana::conv::kPhypDeltaCP,  0);
  skexp.Generate(ph);
  noosc_banff = nsk_integration(noosc_banff, &skexp, s);

  //undo banff
  skexp.GetBANFF()->RestoreNominalMC();
  skexp.Generate(ph);
  noosc_nobanff = nsk_integration(noosc_nobanff, &skexp, s);


  cout << "\n\n\nWith BANFF tuning" << endl;
  nsk_print(osc_banff,input_dm23sq);

  cout << "\n\n\nWithout BANFF tuning" << endl;
  nsk_print(osc_nobanff,input_dm23sq);

  cout << "\n\n\nWith BANFF tuning or osc" << endl;
  nsk_print(noosc_banff,0);

  cout << "\n\n\nWithout BANFF tuning or osc" << endl;
  nsk_print(noosc_nobanff,0);

  //summary table, w & w/o banff & osc
  stringstream streamer;
  string endcolumn  = " & ";
  string endrow     = " \\\\ ";
  int ncol1   = 30;
  int ncolnev = 15;
  int nprecisionnev = 3;
  streamer << "\\begin{table}"   << endl;
  streamer << "\\footnotesize" << endl;
  //streamer << "\\centering"     << endl;
  streamer << "\\begin{tabular}{| l | c c | c c |}" << endl;
  streamer << "\\hline" << endl;
  streamer << " &  "
	   << "\\multicolumn {2}{c|} { {\\bf Osc ($\\Delta{m}_{23}^2$ = " << input_dm23sq << ") } } &  "
	   << "\\multicolumn {2}{c|} { {\\bf No osc} } \\\\"
	   << endl;
  streamer << " &  "
	   << "{\\bf No BANFF}   &  "
	   << "{\\bf BANFF tuned} &  "
	   << "{\\bf No BANFF}   &  "
	   << "{\\bf BANFF tuned}  \\\\"
	   << endl;
  streamer << "\\hline" << endl;

  streamer << fixed;

  streamer << setw(ncol1)  << setfill(' ') << "Total"
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_nobanff.nsk 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_banff.nsk
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_nobanff.nsk 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_banff.nsk
	   << endrow
	   << endl;

  streamer << "\\hline" << endl;
  cout << streamer.str() << endl;
  streamer << setw(ncol1)  << setfill(' ') << "$\\nu_\\mu$ CC"
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_nobanff.nskmucc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_banff.nskmucc
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_nobanff.nskmucc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_banff.nskmucc
	   << endrow
	   << endl;
  streamer << setw(ncol1)  << setfill(' ') << "$\\bar{\\nu_\\mu}$ CC"
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_nobanff.nskmubarcc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_banff.nskmubarcc
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_nobanff.nskmubarcc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_banff.nskmubarcc
	   << endrow
	   << endl;
  streamer << setw(ncol1)  << setfill(' ') << "$\\nu_e$ CC"
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_nobanff.nskecc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_banff.nskecc
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_nobanff.nskecc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_banff.nskecc
	   << endrow
	   << endl;
  streamer << setw(ncol1)  << setfill(' ') << "$\\bar{\\nu_e}$ CC"
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_nobanff.nskebarcc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_banff.nskebarcc
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_nobanff.nskebarcc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_banff.nskebarcc
	   << endrow
	   << endl;
  streamer << setw(ncol1)  << setfill(' ') << "Signal $\\nu_e$ CC"
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_nobanff.nsksigecc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_banff.nsksigecc
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_nobanff.nsksigecc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_banff.nsksigecc
	   << endrow
	   << endl;
  streamer << setw(ncol1)  << setfill(' ') << "NC"
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_nobanff.nsknc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << osc_banff.nsknc
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_nobanff.nsknc 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << noosc_banff.nsknc
	   << endrow
	   << endl;

  streamer << "\\hline" << endl;
  streamer << "\\end{tabular}" << endl;
  streamer << "\\end{table}" << endl;

  cout << "\n\n\n\n\n\n" << streamer.str() << endl;
}


void nsk_print(nsk_t n, double input_dm23sq) {
  //output in form of table 16, TN52
  cout << endl << " dm23sq is " << input_dm23sq << endl;
  cout << "     " << n.nskmucc + n.nskmubarcc + n.nskecc + n.nskebarcc + n.nsksigecc + n.nsknc << " beam total" << endl;
  cout << "     " << n.nskmucc + n.nskmubarcc << " numu+numubar CC" << endl;
  cout << "     " << n.nskecc + n.nskebarcc << " nue+nuebar beam CC" << endl;
  cout << "     " << n.nsknc << " all NC" << endl;
  cout << "     " << n.nsksigecc << " nue signal CC" << endl;

  //output a bit more broken down
  cout << endl << " dm23sq is " << input_dm23sq << endl;
  cout << "     " << n.nskmucc + n.nskmubarcc + n.nskecc + n.nskebarcc + n.nsksigecc + n.nsknc << " beam total" << endl;
  cout << "     " << n.nskmucc << " numu CC" << endl;
  cout << "     " << n.nskmubarcc << " numubar CC" << endl;
  cout << "     " << n.nskecc << " nue beam CC" << endl;
  cout << "     " << n.nskebarcc << " nuebar beam CC" << endl;
  cout << "     " << n.nsknc << " all NC" << endl;
  cout << "     " << n.nsksigecc << " nue signal CC" << endl;

  //the total
  cout << endl << "Number of predicted events at SK is: " << n.nsk << "     = " 
       << n.nskmucc + n.nskmubarcc + n.nskecc + n.nskebarcc + n.nsksigecc + n.nsknc << "     = " 
       << n.nskmc << endl;
}

nsk_t nsk_integration(nsk_t n, SuperkExpectation *skexp, conv::SuperkSample_t s) {
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  // integrate Ereco distribution
  double n.nsk = skexp->CurrPredictedErecoSpectrum(s)->Integral();
  double n.nskmc = 0, n.nskmucc = 0, n.nskmubarcc = 0, n.nskecc = 0, n.nskebarcc = 0, n.nsksigecc = 0, n.nsknc = 0;
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined)
      continue;
    n.nskmc += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();

    //split into categories (all NC, CC numu, CC numubar, CC nue, CC nuebar, CC signue)
    if(util::IsNC(m))
      n.nsknc += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
    else {
      if(util::IsNuMu(m))
	n.nskmucc += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
      else if(util::IsNuMuBar(m))
	n.nskmubarcc += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
      else if(util::IsNuE(m))
	n.nskecc += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
      else if(util::IsNuEBar(m))
	n.nskebarcc += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
      else if(util::IsSignalNuE(m))
	n.nsksigecc += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
    }//CC?
  }
  return n;
}
