/*

root2k -l -b -q '$T2KOSCANA/macros/misc/automated_caption_maker.C("./images/spectrum_tweaked/unoscillated/","ereco_twk_","tweakspectra",0,0,0,0,0,0)'
root2k -l -b -q '$T2KOSCANA/macros/misc/automated_caption_maker.C("./images/spectrum_tweaked/oscillated/","ereco_twk_","tweakspectra",-99,-99,-99,-99,-99,-99)'

*/

void automated_caption_maker(
  const char * graphicslocation,                                  /* folder where the plots are */
  const char * graphicsprefix,                                    /* prefix on the plots (end is assumed to be AsStringCompact) */
  const char * tag              = "tweakspectra",                 /* output tex file name prefix */
  double       input_mixing12   = -99,                            /* */
  double       input_mixing13   = -99,                            /* */
  double       input_mixing23   = -99,                            /* */
  double       input_deltacp    = -99,                            /* */
  double       input_dm12sq     = -99,                            /* */
  double       input_dm23sq     = -99                             /* */
)
{

  using namespace std;
  using namespace t2kana;

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


  //loop over systematic parameters
  for(conv::PhysParamId_t p = conv::kNuisSuperkEnergyScale; p <= conv::kNuisFSI005; p++) {

    if(p == conv::kNuisXSecENuShape || p == conv::kNuisXSecNormNCcoh)
      continue;

    string syscompact = util::AsStringCompact(p);
    string syspaper   = util::AsStringPaperTex(p);

    //front matter
    streamer << "\\begin{figure}" << endl;
    streamer << "\\begin{center}" << endl;
    streamer << "\\resizebox{14cm}{!}{" 
	     << Form("\\includegraphics{%s%s%s}}",graphicslocation,graphicsprefix,syscompact.c_str()s) << endl;

    //caption
    streamer << "\\caption" << endl;
    streamer << "{ Effect of " << syspaper.c_str() << " for ";
    string oscstr = "";
    if(input_mixing12 < conv::kASmallNum && input_mixing13 < conv::kASmallNum 
       && input_mixing23 < conv::kASmallNum && input_dm12sq < conv::kASmallNum
       && input_dm23sq < conv::kASmallNum && input_dm12sq < conv::kASmallNum) {
      streamer << "no oscillations (all mixing angles set to 0)";
      oscstr = "noosc";
    }
    else {
      streamer << "for 3 flavour oscillations,"
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
    streamer << Form("\\label{fig:%s_%s%s}",oscstr.c_str(),graphicsprefix,syscompact.c_str()) << endl;


    //back matter
    streamer << "\\end{center}" << endl;
    streamer << "\\end{figure}" << endl;
    streamer << "\\clearpage" << endl << endl;

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

