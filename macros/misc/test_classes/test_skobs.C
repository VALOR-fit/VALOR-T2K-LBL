//
// ------------------------------------------------------------------------------------------------
// Run all the SKObs methods to test that they work correctly
//
// Tom Dealtry
//
// ------------------------------------------------------------------------------------------------
//

void test_skobs
(
 bool test_read = false,
 bool test_generate = true,
 bool test_num_of_bins = false,
 bool test_num_of_events = false,
 const char * tag                     = "testSKOBS",    /* tag for output file */
 int          run                     = 0,                /* run number */
 int          first_experiment        = 0,                       /* */
 int          last_experiment         = 0,                       /* */
 const char * xmlfile = "test/twod_toys.xml",
 bool         using_xml_file          = true,                                /* if false, skobs reads xmlfile as a root file */
 const char * config_file      = "analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml"
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  if(test_read == test_generate) {
    cerr << "Exactly one of test_read and test_generate must be turned on" << endl;
    exit(-1);
  }

  // Analysis configuration
  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file));
  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile));

  FileStat_t dummy;

  // Check for config file, and XML data file.
  // If not found in inputs or config, see if it works as a normal path
  if (gSystem->GetPathInfo(config.Data(), dummy)) 
  {
    config = config_file;
    if (gSystem->GetPathInfo(config.Data(), dummy)) {
      cerr<<"Config file not found: "<<config_file<<endl;
      exit(-1);
    }
  }
  
  if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) 
  {
    xmlfile2 = xmlfile;
    if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
      cerr<<"Data input file not found: "<<xmlfile<<endl;
      exit(-1);
    }
  }

  int nexperiments = last_experiment - first_experiment + 1; 
  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
  cout << "   - Reading experiments from: " << xmlfile2 << endl; 
  cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl; 

  TString outfilename = Form("%s-run_%d_%d_%d.root", tag, run, first_experiment, last_experiment);
  TFile fout(outfilename.Data(), "RECREATE");
  fout.Close();

  // Loop over toy-MC experiments
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {
    cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

    // get toy MC with/without systematics
    
    SuperkObservation * skobs = new SuperkObservation;
    SuperkExpectation * skexp;

    if(test_read)
      skobs->Read(xmlfile2.Data(),iexperiment,using_xml_file);

    if(test_generate) {
      skexp = new SuperkExpectation;
      skexp->Read(config.Data());
      PhysHypothesis ph(kHypOsc3F);
      ph.SetCurrValue(kPhypMixing12, kT2KBestFitSinsqTheta12);
      ph.SetCurrValue(kPhypMixing23, kT2KBestFitSinsqTheta23_NH);
      ph.SetCurrValue(kPhypMixing13, kT2KBestFitSinsqTheta13_NH);
      ph.SetCurrValue(kPhypDm12sq,   kT2KBestFitDmsq12);
      ph.SetCurrValue(kPhypDm23sq,   kT2KBestFitDmsq23_NH);
      ph.SetCurrValue(kPhypCPV,      kT2KBestFitDeltaCP_NH);

      bool fluctuate = false;
      skobs->Generate(skexp, ph, fluctuate);

      double chisq = skobs->Statistic(skexp);
      cout << "Total chisq = " << chisq << endl;
    }//test_generate

    if(test_num_of_bins)
      cout << "NBins = " << skobs->NumOfBins() << endl;

    if(test_num_of_events) {
      for (int is = 0 ; is < conv::kns ;  is ++ ) {
	SuperkSample_t s = kSKSamples[is];
	if(s == kSKsUndefined) continue;
	cout << "NEv = " << skobs->NumOfEvents(s) << "\t" << util::AsStringCompact(s) << endl;
	skobs->NumOfEvents(-1,-1,s);
	for(int ir = 0; ir < util::NumErecoBins(s); ir++)
	  skobs->NumOfEvents(util::ErecoBinEdges(s,ir) + conv::kASmallNum, -1, s);
	for(int iq = 0; iq < util::NumQ2Bins(s); iq++)
	  skobs->NumOfEvents(-1, util::Q2BinEdges(s,iq) + conv::kASmallNum, s);
	for(int ir = 0; ir < util::NumErecoBins(s); ir++)
	  for(int iq = 0; iq < util::NumQ2Bins(s); iq++)
	    skobs->NumOfEvents(util::ErecoBinEdges(s,ir) + conv::kASmallNum,
			       util::Q2BinEdges(s,iq) + conv::kASmallNum, s);
      }//is
    }//test_num_of_events

    cout << "________________________ READ, NOW LETS WRITE ____________" << endl;

    skobs->Write(outfilename.Data(), "obs", "hist", "update");
    if(test_generate)
      skexp->Write(outfilename.Data(), "exp", "hist", "update");

    // Clean-up
    delete skobs;
    delete skexp;
  }//iexperiment
}
