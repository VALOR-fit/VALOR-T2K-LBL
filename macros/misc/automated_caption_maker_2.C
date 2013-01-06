/*
used to make the captions for the appendix with the star plots
*/

void automated_caption_maker_2(
  const char * graphicslocation = "./images/star_plot/",          /* folder where the plots are */
  const char * tag              = "star_plots",                   /* output tex file name prefix */
  double       input_mixing12   = -99,                            /* */
  double       input_mixing13   = -99,                            /* */
  double       input_deltacp    = -99,                            /* */
  double       input_dm12sq     = -99                             /* */
)
{

  using namespace std;
  using namespace t2kana;

  // Set inputs
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;


  string output_file (Form("%s.tex",tag));

  stringstream streamer;
  streamer << "\\clearpage" << endl << endl;

  //string of systematics for filenames & labels
  vector<string> systematics;
  systematics.push_back("important");
  systematics.push_back("sk");
  systematics.push_back("banffnumu");
  systematics.push_back("banffnumubar");
  systematics.push_back("banffxsec");
  systematics.push_back("xsec");
  systematics.push_back("fsi");
  systematics.push_back("osc");

  //string of systematic explanations, for caption
  vector<string> sysexplanation;
  sysexplanation.push_back("parameters used in the 2+8 fit");
  sysexplanation.push_back("SuperK detector parameters");
  sysexplanation.push_back("BANFF $\\nu_\\mu$ flux parameters");
  sysexplanation.push_back("BANFF $\\bar{\\nu_\\mu}$ flux parameters");
  sysexplanation.push_back("BANFF cross section parameters");
  sysexplanation.push_back("uncorrelated cross section parameters");
  sysexplanation.push_back("FSI+SI parameters");
  sysexplanation.push_back("non-23 sector oscillation parameters");

  vector<double> mixing23;
  mixing23.push_back(0.8);
  mixing23.push_back(0.8);
  mixing23.push_back(0.8);
  mixing23.push_back(0.9);
  mixing23.push_back(0.9);
  mixing23.push_back(0.9);
  mixing23.push_back(0.95);
  mixing23.push_back(0.95);
  mixing23.push_back(0.95);
  mixing23.push_back(1.0);
  mixing23.push_back(1.0);
  mixing23.push_back(1.0);
  mixing23.push_back(1.0);
  mixing23.push_back(1.0);

  vector<double> dm23sq;
  dm23sq.push_back(0.002);
  dm23sq.push_back(0.0024);
  dm23sq.push_back(0.003);
  dm23sq.push_back(0.0022);
  dm23sq.push_back(0.0024);
  dm23sq.push_back(0.0026);
  dm23sq.push_back(0.0022);
  dm23sq.push_back(0.0024);
  dm23sq.push_back(0.0026);
  dm23sq.push_back(0.002);
  dm23sq.push_back(0.0022);
  dm23sq.push_back(0.0024);
  dm23sq.push_back(0.0026);
  dm23sq.push_back(0.003);

  for(int isys = 0; isys < systematics.size(); isys++) {
    string systematic = systematics[isys];
    string sysexplain = sysexplanation[isys];
    string filesuffix = "_testzoom";
    if(systematic == "important")
      filesuffix = "_testbig";

    string sysexplainupper = sysexplanation[isys];
    sysexplainupper[0] = toupper(sysexplain[0]);
    streamer << endl << "\\subsection{" << sysexplainupper << "}" << endl;

    for(int iosc = 0; iosc < mixing23.size(); iosc++) {

      double input_mixing23 = mixing23[iosc];
      double input_dm23sq = dm23sq[iosc];

      //front matter
      if(iosc % 4 == 0) {
	streamer << "\\begin{figure*}[htb]" << endl;
      }
      streamer << "\\begin{minipage}[t]{0.48\\linewidth}" << endl;
      streamer << "\\center" << endl;

      //image
      streamer << Form("\\includegraphics[width=18pc]{{%s%s-%.6f-%.6f%s}.eps}",graphicslocation,systematic.c_str(),input_mixing23,input_dm23sq,filesuffix.c_str()) << endl;

      //caption
      streamer << "\\caption" << endl;
      streamer << "{ Effect of " << sysexplain.c_str() << " on the best fit point, using a 2+0 fit, with a single systematic changed by $\\pm$1 sigma. "
	       << "The ellipse shown is one quarter the size of the statistical error. For ";
      string oscstr = "";
      if(input_mixing12 < conv::kASmallNum && input_mixing13 < conv::kASmallNum 
	 && input_mixing23 < conv::kASmallNum && input_dm12sq < conv::kASmallNum
	 && input_dm23sq < conv::kASmallNum && input_dm12sq < conv::kASmallNum) {
	streamer << "no oscillations (all mixing angles set to 0)";
	oscstr = "noosc";
      }
      else {
	streamer << "3 flavour oscillations,"
		 << " using oscillation parameters: " 
		 << " $\\sin^22\\theta_{12} = $"  << input_mixing12
		 << "; $\\sin^22\\theta_{13} = $" << input_mixing13
		 << "; $\\sin^22\\theta_{23} = $" << input_mixing23
		 << "; $\\delta_{CP} = $        " << input_deltacp
		 << "; $\\Delta m^2_{12} = $    " << input_dm12sq
		 << "; $\\Delta m^2_{32} = $    " << input_dm23sq;
	oscstr = "osc";
      }
      streamer << ". }" << endl;

      //label
      streamer << Form("\\label{fig:apdx_star_%s_%.2f_%.4f}",systematic.c_str(),input_mixing23,input_dm23sq) << endl;

      //back matter
      streamer << "\\end{minipage}" << endl;
      if(iosc % 4 == 3 || iosc == mixing23.size()-1) {
	streamer << "\\end{figure*}" << endl;
	streamer << "\\clearpage" << endl;
	//if(iosc %4 == 1 || iosc == mixing23.size() -1) {
      }
      else if(iosc %4 == 0 || iosc %4 == 2){
	streamer << "\\hspace{0.06\\linewidth}" << endl;
      }
      else if(iosc %4 == 1) {
	streamer << "\\vfill" << endl;
      }


    }//isys
  }//iosc

  cout << streamer.str() << endl;

  // write the table
  long size;
  // get size of file
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);
  ofstream outfile (output_file.c_str(),ios_base::out);
  // write to outfile
  outfile.write (streamer.str().c_str(),size);
  outfile.close();
}

