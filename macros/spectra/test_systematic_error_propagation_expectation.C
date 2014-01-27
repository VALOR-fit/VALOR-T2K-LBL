#include <TXMLEngine.h>
#include <iomanip.h>

test_systematic_error_propagation_expectation(
  bool         use_tex      = true,
  const char * output_file  = "systematic_tables_expectation.tex",
  const char * config_file  = "analysis_configuration.xml",    /* */
  bool         for_tn       = false
)
{
  using namespace std;
  using namespace t2kana;
  
  // Parse the input XML configuration file
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config << endl;
    delete xml;
    exit(1);
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

  // Get the expected systematic variations

  //SK
  double SKDetEnergyScale1RmuRelErr = util::XML2Dbl(xml,root_node,
     "detector_systematics/superk/energy_scale_relative_error/one_mu_ring");

  //BANFF
  TMatrixD *BANFFCovMtrx = util::XML2TMatrixD(xml,root_node, "banff/covariance_matrix");
  TVectorD *BANNFBestFit = util::XML2TVectorD(xml,root_node, "banff/best_fit_values");

  //SKDetFSI
  TMatrixD *SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node, "skdetfsi/covariance_matrix");

  //XSec
  double NormCCcohRelErr = util::XML2Dbl(xml,root_node,"xsec/normalisations/cccoh");
  double NormNC1piCRelErr = util::XML2Dbl(xml,root_node,"xsec/normalisations/nc1piC");
  double NormNCcohRelErr = util::XML2Dbl(xml,root_node,"xsec/normalisations/nccoh");
  double NormNCothRelErr = util::XML2Dbl(xml,root_node,"xsec/normalisations/ncoth");
  double NormNuEToNuMuRelErr = util::XML2Dbl(xml,root_node,"xsec/normalisations/numu_to_nue");
  double NormNuToAntiNuRelErr = util::XML2Dbl(xml,root_node,"xsec/normalisations/nu_to_antinu");

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  //print out
  string pctsymbol  = (use_tex) ? " \\% "   : "% ";
  string endcolumn  = (use_tex) ? " & "     : "" ;
  string endrow     = (use_tex) ? " \\\\ "  : "" ;
  int nprecisionnev = 4; // <-- for number of events
  int nprecisionpct = 2; // <-- for % changes
  int ncol1   = 80;
  int ncolnev = 15;
  int ncolfrc =  8;

  stringstream streamer;

  if(use_tex) {
    streamer << "\\documentclass[10pt,a4paper]{article}\n" << endl;
    streamer << "\\usepackage{amsmath}\n" << endl;
    streamer << "\\usepackage{graphicx}\n" << endl;
    streamer << "\\usepackage{dcolumn}\n" << endl;
    streamer << "\\usepackage[cm]{fullpage}\n" << endl;
    streamer << "\\begin{document}\n" << endl;

    streamer << "\\begin{table}"   << endl;
    streamer << "\\scriptsize" << endl;
    streamer << "\\centering"     << endl;
    streamer << "\\begin{tabular}{| l | r |}" << endl;
    streamer << "\\hline" << endl;
    streamer << setw(ncol1)  << setfill(' ') << "{\\bf Systematic}  " << endcolumn
	     << setw(ncolnev)  << setfill(' ') << "{\\bf Expected change} " << endrow
	     << endl;
    streamer << "\\hline" << endl;
  }

  string systematic;
  systematic = "$f_{E; r}^{SK}$";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << 0 << pctsymbol << " (movement in Ereco)"
	   << endrow << endl;
  for (int i = 0; i < conv::kNumOfParamsInBANFFfit; i++) {
    //if(i == 16 || i == 17) continue; //skip MAQE, MARES
    systematic = Form("$f_{%d; t, r}^{banff}$ ",i);
    if(!for_tn)
      systematic += util::AsStringTex(conv::kNuisBANFF00+(conv::PhysParamId_t)i);
    streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	     << sqrt(BANFFCovMtrx(i,i)) * 100 / BANNFBestFit(i) << pctsymbol;
    streamer << endrow << endl;
  }
  systematic = "$f_{CCcoh; t}$";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormCCcohRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "$f_{NC1\\pi^\\pm; t}$";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormNC1piCRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "$f_{NCcoh; t}$";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormNCcohRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "$f_{NCoth; t}$";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormNCothRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "$f_{CC\\nu_e; t}$";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormNuEToNuMuRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "$f_{CC\\bar{\\nu}; t}$";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormNuToAntiNuRelErr * 100 << pctsymbol
	   << endrow << endl; NormNuToAntiNuRelErr;
  for (int i = 0; i < conv::kNumOfSKDetFSIParams; i++) {
    systematic = Form("$f_{%d; t, r}^{SKFSI}$ ",i);
    if(!for_tn)
      systematic += util::AsStringTex(conv::kNuisSKDetFSI000+(conv::PhysParamId_t)i);
    streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	     << sqrt(SKDetFSICovMtrx(i,i)) * 100 << pctsymbol;
    streamer << endrow << endl;
  }
  if(!use_tex) {
     streamer << setw(150)  << setfill('-') << " " << endl;
  }

  if(use_tex) {
     streamer << "\\hline" << endl;
     streamer << "\\end{tabular}" << endl;
     streamer << "\\caption{Summary table, showing the magnitude of errors on the systematic parameters. The numbers should be compared to the changes in $N_{SK}$ in the following tables. Systematics with complicated energy dependence, treated with response functions, are not shown.}" << endl;
     streamer << "\\label{tab:expected}" << endl;
     streamer << "\\end{table}" << endl;
     streamer << "\\clearpage" << endl;
     streamer << "\\end{document}" << endl;
  }

  cout << streamer.str() << endl;

  // write the table?
  if(output_file != "") {
    long size;

    // get size of file
    streamer.seekg(0,ifstream::end);
    size=streamer.tellg();
    streamer.seekg(0);

    ofstream outfile (output_file,ios_base::app);

    // write to outfile
    outfile.write (streamer.str().c_str(),size);
  
    outfile.close();
  }

}
