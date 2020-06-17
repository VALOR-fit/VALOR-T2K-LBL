
using namespace valor;
using namespace valor::conv;
using namespace valor::util;
using namespace std;

TSpline3 * RandomSpline (double low_knot, double high_knot, double sigma)
{
  static double x[2] = {low_knot,high_knot};
  static double y[2];

  int slope = gRandom->Integer(2)*2 -1;
  double low_diff  = gRandom->Uniform(-3 * sigma,3*sigma);
  double high_diff = gRandom->Uniform(-3 * sigma,3*sigma);

  y[0] = 1 - (low_diff  * slope); 
  y[1] = 1 + (high_diff * slope); 

  return new TSpline3("FAKESPLINE",x,y,2);
}


void make_fake_splines(const char * version_number = "2016v0", const char * config = "config.xml", const char * detector = "dune_nd")
{

  int n_splines = 0;
  TFile f(Form("fake_splines_%s.root",version_number),"RECREATE");

  int izone  = 0;
  int itrue3 = 0;
  int ireco3 = 0;

  // Read config to populate the ParameterPool
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config);

  std::vector<conv::PhysParamId_t> params_with_splines;
  for (const Parameter & param_search : ParameterPool::Instance()->GetParams()) {
    if(param_search.GetType() == conv::kSystResponse) {
      params_with_splines.push_back(param_search.GetID());
    }
  }

  for(int ibeam = 0; ibeam < kNBeamConf; ibeam++) {
    BeamConf_t b = kBeamConf[ibeam];
    assert(b != kbUndefined);
    for(int isample = 0; isample < kNSamples; isample++) {
      Sample_t s = kSample[ibeam][isample];
      assert(s != ksUndefined);
      
      cout<<util::AsString(b)<<"\t"<<util::AsString(s)<<endl;

      RFL rfl(b,s,params_with_splines);

      const int n_enu_true_bins = util::NumEtrueBinsXsec (s);
      const int n_enu_reco_bins = util::NumErecoBinsXsec (s, 0);
      const int n_y_reco_bins   = util::NumYrecoBinsXsec (s, 0);
        
      for (const PhysParamId_t p : params_with_splines) {
        assert(p != kPhypUndefined);

        for(int imode = 0; imode < kNReacModes; imode++) {
          ReacMode_t m = kReactionModes[isample][imode];
          if(m == kmUndefined) continue;

          const Parameter & param = ParameterPool::Instance()->GetParam(p);
          conv::Detector_t d = util::DetectorFromString(detector);
          if (! util::IsAppliedTo(param.GetApplyMask(),b,s,d,m,conv::kntInclusive) ) {continue;}

          const int n_q2_true_bins = util::NumQ2trueBinsXsec (s, m);
          
          KinematicBin kb;
          kb.zone = 0;
          kb.true3 = 0;
          kb.reco3 = 0;

          for (int itrue1 = 0 ; itrue1 < n_enu_true_bins ; itrue1++) {
            kb.true1 = itrue1;
            for (int itrue2 = 0 ; itrue2 < n_q2_true_bins ; itrue2++) {
              kb.true2 = itrue2;
              for (int ireco1 = 0 ; ireco1 < n_enu_reco_bins ; ireco1++) {
                kb.reco1 = ireco1;
                for (int ireco2 = 0 ; ireco2 < n_y_reco_bins ; ireco2++) {
                  kb.reco2 = ireco2;

                  TString spl_name = rfl.BuildSplineName(m,p,kb);

                  TSpline3 * spl = RandomSpline(0,2,0.3);
                  spl->SetName(spl_name);
                  spl->Write();
                  delete spl;
                  
                  n_splines++;
                  if ((n_splines%1000)==0) {
                    cout<<n_splines<<endl;
                  }
                } //reco2
              } //reco1
            }//true2
          }// true1
        }// mode
      }// sample
    } //beam
  }// parameter
  cout<<"Total: "<<n_splines<<endl;
  f.Close();
  return;
}
