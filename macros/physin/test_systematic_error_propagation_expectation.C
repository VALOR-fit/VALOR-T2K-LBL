#include <TXMLEngine.h>
#include <iomanip.h>

test_systematic_error_propagation_expectation(
  bool         use_tex      = false,
  const char * output_file  = "",
  const char * config_file  = "analysis_configuration.xml"    /* */
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
  double SKDetEnergyScale1RmuRelErr = 
    util::XML2Dbl(xml,root_node,
     "detector_systematics/superk/energy_scale_relative_error/one_mu_ring");
  double SKDetEff1RmuCCQEOthRelErr =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/other_error");
  TMatrixD *SKDetEff1RmuCCQERCCovarRelErr = 
    util::XML2TMatrixD(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/rc_error");
  double SKDetEff1RmuCCnonQERelErr   =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccnonqe");
  double SKDetEff1RmuNCRelErr        =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/nc");
  double SKDetEff1RmuNuECCRelErr =
    util::XML2Dbl(xml,root_node, 
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/nuecc");

  //BANFF
  TMatrixD *BANFFCovMtrx = 
    util::XML2TMatrixD(xml,root_node,
      "banff/covariance_matrix");

  //XSec
  double NormCCQE1RelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/ccqee1");
  double NormCCQE2RelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/ccqee2");
  double NormCC1pi1RelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/cc1pie1");
  double NormCCcohRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/cccoh");
  double NormNCcohRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/nccoh");
  double NormNCothRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/ncoth");
  double NormNuEToNuMuRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/numu_to_nue");

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
     streamer << "\\begin{table}"   << endl;
     streamer << "\\scriptsize" << endl;
     //streamer << "\\centering"     << endl;
     streamer << "\\begin{tabular}{| l | r |}" << endl;
     streamer << "\\hline" << endl;
     streamer << setw(ncol1)  << setfill(' ') << "{\\bf Systematic}  " << endcolumn
	      << setw(ncolnev)  << setfill(' ') << "{\\bf Expected change} " << endrow
	      << endl;
     streamer << "\\hline" << endl;
  }
 
  string systematic;

  for (int i = 0; i < 3; i++) {
    systematic = Form("sk 1Rmu ccqe ring counting eff %d",i);
    streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	     << TMath::Sqrt(SKDetEff1RmuCCQERCCovarRelErr(i,i)) * 100 << pctsymbol
	     << endrow << endl;
  }
  systematic = "sk 1Rmu ccqe other eff";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << SKDetEff1RmuCCQEOthRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "sk 1Rmu ccnonqe eff";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << SKDetEff1RmuCCnonQERelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "sk 1Rmu nc eff";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << SKDetEff1RmuNCRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "sk 1Rmu nuecc eff";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << SKDetEff1RmuNuECCRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "sk energy scale";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << 0 << pctsymbol << " (movement in Ereco)"
	   << endrow << endl;
  for (int i = 0; i < conv::kNumOfParamsInBANFFfit; i++) {
    systematic = Form("banff %d",i);
    streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	     << BANFFCovMtrx(i,i) * 100 << pctsymbol;
    streamer << endrow << endl;
  }
  systematic = "Norm CCQE E1";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormCCQE1RelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "Norm CCQE E2";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormCCQE2RelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "Norm CC1pi E1";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormCC1pi1RelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "Norm CC coh";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormCCcohRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "Norm NC coh";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormNCcohRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "Norm NC oth";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormNCothRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "Norm NuE to NuMu";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << NormNuEToNuMuRelErr * 100 << pctsymbol
	   << endrow << endl;

  if(!use_tex) {
     streamer << setw(150)  << setfill('-') << " " << endl;
  }

  if(use_tex) {
     streamer << "\\hline" << endl;
     streamer << "\\end{tabular}" << endl;
     streamer << "\\label{tab:expected}" << endl;
     streamer << "\\end{table}" << endl;
     streamer << "\\clearpage" << endl;
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
