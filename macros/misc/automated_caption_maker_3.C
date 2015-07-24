/*
used to create the captions for the etrue binning appendix subsection


root2k -l -b -q '$T2KOSCANA/macros/misc/automated_caption_maker_3.C("./images/binning/oscillation/")'

*/

void automated_caption_maker_3(
  const char * graphicslocation,                                  /* folder where the plots are */
  const char * graphicsprefix   = "plot_erec_int_",                                    /* prefix on the plots (end is assumed to be AsStringCompact) */
  const char * tag              = "etruebinning",                 /* output tex file name prefix */
  double       input_mixing12   = -99,                            /* */
  double       input_mixing13   = -99,                            /* */
  double       input_mixing23   = -99,                            /* */
  double       input_deltacp    = -99,                            /* */
  double       input_dm12sq     = -99,                            /* */
  double       input_dm23sq     = -99                             /* */
)
{

  using namespace std;
  using namespace valor;

  // Set inputs
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

  string oscparams = util::AsString(conv::kHypOsc3F,
	  input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  string output_file (Form("%s-%s.tex",tag,oscparams.c_str()));

  stringstream streamer;

  iplot = 0;
  //loop over systematic parameters
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[0][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::IsNC(m)) continue;

    string modecompact = util::AsStringCompact(m);
    string modecompactJ = util::AsStringCompact_Japan(m);
    string modelong = util::AsStringTex(m);

    //front matter
    if(iplot % 4 == 0) {
      streamer << "\\begin{figure*}[htb]" << endl;
    }
    streamer << "\\begin{minipage}[t]{0.48\\linewidth}" << endl;
    streamer << "\\center" << endl;

    //image
    streamer << Form("\\includegraphics[width=18pc]{%s%s%s_fine}",graphicslocation,graphicsprefix,modecompactJ.c_str()) << endl;

    //caption
    streamer << "\\caption" << endl;

    streamer << " {Ratio of SuperK MC $E_{reco}$ spectra, "
	     << "for mode " << modelong.c_str()
	     << " made with 84 bin $E_{true}$ and 780 bin $E_{true}$, for an exposure of 3.010$\times$10$^{20}$ POT. "
	     << "The number of events was calculated using the BANFF-tuned MC templates.";

    streamer << "Using oscillation parameters: " 
	     << " $\\sin^22\\theta_{12} = $"  << input_mixing12
	     << "; $\\sin^22\\theta_{13} = $" << input_mixing13
	     << "; $\\sin^22\\theta_{23} = $" << input_mixing23
	     << "; $\\delta_{CP} = $        " << input_deltacp
	     << "; $\\Delta m^2_{21} = $    " << input_dm12sq << " eV^2/c^4"
	     << "; $|\\Delta m^2_{32}| = $    " << input_dm23sq << " eV^2/c^4";
    streamer << ". }" << endl;

    //label
    streamer << Form("\\label{fig:etruebin:%s}",modecompact.c_str()) << endl;


    //back matter
    streamer << "\\end{minipage}" << endl;
    if(iplot % 4 == 3){// || iplot == mixing23.size()-1) {
      streamer << "\\end{figure*}" << endl << endl;
      streamer << "\\clearpage" << endl;
      //if(iplot %4 == 1 || iplot == mixing23.size() -1) {
    }
    else if(iplot %4 == 0 || iplot %4 == 2){
      streamer << "\\hspace{0.06\\linewidth}" << endl;
    }
    else if(iplot %4 == 1) {
      streamer << "\\vfill" << endl;
    }
    iplot++;

  }

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

