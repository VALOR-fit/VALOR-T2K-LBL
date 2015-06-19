/*
This macro makes the plots of P(numu->nue) against P(numubar->nuebar) for T2K and NOvA and P(numu->nue) for T2K against 
P(numu->nue) for NOvA that were shown at the T2K collaboration meeting on 25th April 2012. These slides are at

http://www.t2k.org/meet/collab/201204/talks/wedam/DeltaCP/view

If useMatterProb is true, the plots are made including matter effects as in slides 3, 6, 9-10 and 19. 

If useMatterProb is false, the plots are made using vacuum probabilities as in slides 5 and 8. The vacuum version of 
slide 19 is also made, but the plots with errors on slides 9-10 are not made. 

Nick Grant
*/

void FormatGraphLong(TGraph *graph, int colour, int width, int style);

void FormatGraph(TGraph *graph, int colour, int width, int style);

void plot_deltacp_hierarchy_comp_t2k_nova(
bool useMatterProb = true,
double sinsq2theta12 = 0.87,
double sinsq2theta13 = 0.1,
double sinsq2theta23 = 1.0,
double dmsq12 = 7.6e-5,
double modDmsq32 = 2.32e-3,
double sinsq2theta13Error = 0.02, //current value 0.02, Daya Bay + Reno 0.01, T2K LOI 0.006
//Don't need a check on sinsq2theta23 becoming > 1.0 with the addition of its error 
//since NuOsc3Flavour restricts it to the physically allowed region between 0 and 1
double sinsq2theta23PlusError = 0.0,  //current value 0.0 (since best fit is 1.0), T2K LOI 0.01
double sinsq2theta23MinusError = 0.06,  //current value 0.06, T2K LOI 0.01 
double modDmsq32PlusError = 0.12e-3, //current value 0.12e-3, T2K LOI 0.1e-3
double modDmsq32MinusError = 0.08e-3,  //current value 0.08e-3, T2K LOI 0.1e-3 
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  // Global options
  Opt::Instance()->gOscMixingParamChoice = kNuMixSinsq2Theta;

  NuOsc3Flavour * nuOsc = new NuOsc3Flavour;

  const int npoints = 5000;

  double energy, deltaCP, oscProb;
  const double pi = 3.14159265;

  double deltaCPpoints[4] = {0.0, 0.5*pi, pi, 1.5*pi};

  double t2kNumuToNueNH[npoints];
  double t2kNumuToNueIH[npoints];
  double t2kNumubarToNuebarNH[npoints];
  double t2kNumubarToNuebarIH[npoints];
  double novaNumuToNueNH[npoints];
  double novaNumuToNueIH[npoints];
  double novaNumubarToNuebarNH[npoints];
  double novaNumubarToNuebarIH[npoints];

  double t2kOscProbNH[4];
  double t2kOscProbBarNH[4];
  double t2kOscProbIH[4];
  double t2kOscProbBarIH[4];
  double novaOscProbNH[4];
  double novaOscProbBarNH[4];
  double novaOscProbIH[4];
  double novaOscProbBarIH[4];

  double t2kErrorS13NHPlus[4];
  double t2kErrorS13NHMinus[4];
  double t2kErrorS23NHPlus[4];
  double t2kErrorS23NHMinus[4];
  double t2kErrorDmsq32NHPlus[4]; 
  double t2kErrorDmsq32NHMinus[4];
  double t2kErrorS13IHPlus[4];
  double t2kErrorS13IHMinus[4];
  double t2kErrorS23IHPlus[4];
  double t2kErrorS23IHMinus[4];
  double t2kErrorDmsq32IHPlus[4]; 
  double t2kErrorDmsq32IHMinus[4]; 

  double novaErrorS13NHPlus[4];
  double novaErrorS13NHMinus[4];
  double novaErrorS23NHPlus[4];
  double novaErrorS23NHMinus[4];
  double novaErrorDmsq32NHPlus[4]; 
  double novaErrorDmsq32NHMinus[4];
  double novaErrorS13IHPlus[4];
  double novaErrorS13IHMinus[4];
  double novaErrorS23IHPlus[4];
  double novaErrorS23IHMinus[4];
  double novaErrorDmsq32IHPlus[4]; 
  double novaErrorDmsq32IHMinus[4];

  //------------- Calculate P(numu->nue) and P(numubar->nuebar) for different values of deltaCP for T2K ----------------

  energy = 0.6;

  nuOsc->SetMassHierarchy(kNmhNormal);

  for(int i=0; i<npoints; i++)
    {
      deltaCP = i * 6.283185 / npoints;
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32, deltaCP);
      if(useMatterProb)
	{
	  t2kNumuToNueNH[i] = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
	  t2kNumubarToNuebarNH[i] = nuOsc->CalcMatterOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
      else
	{ 
	  t2kNumuToNueNH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMu, kPdgNuE, energy);
	  t2kNumubarToNuebarNH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
    }

  for(int i=0; i<4; i++)
    {
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      if(useMatterProb)
        {
          t2kOscProbNH[i] = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
          t2kOscProbBarNH[i] = nuOsc->CalcMatterOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
        }
      else
        {
          t2kOscProbNH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMu, kPdgNuE, energy);
          t2kOscProbBarNH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
        }

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13+sinsq2theta13Error, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorS13NHPlus[i] = fabs(oscProb - t2kOscProbNH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13-sinsq2theta13Error, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorS13NHMinus[i] = fabs(oscProb - t2kOscProbNH[i]);

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23+sinsq2theta23PlusError, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorS23NHPlus[i] = fabs(oscProb - t2kOscProbNH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23-sinsq2theta23MinusError, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorS23NHMinus[i] = fabs(oscProb - t2kOscProbNH[i]);

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32+modDmsq32PlusError, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorDmsq32NHPlus[i] = fabs(oscProb - t2kOscProbNH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32-modDmsq32MinusError, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorDmsq32NHMinus[i] = fabs(oscProb - t2kOscProbNH[i]);

    }

  nuOsc->SetMassHierarchy(kNmhInverted);

  for(int i=0; i<npoints; i++)
    {
      deltaCP = i * 6.283185 / npoints;
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32, deltaCP);
      if(useMatterProb)
	{ 
	  t2kNumuToNueIH[i] = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
	  t2kNumubarToNuebarIH[i] = nuOsc->CalcMatterOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
      else
	{
	  t2kNumuToNueIH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMu, kPdgNuE, energy);
	  t2kNumubarToNuebarIH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
    }

  for(int i=0; i<4; i++)
    {
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      if(useMatterProb)
	{
	  t2kOscProbIH[i] = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
	  t2kOscProbBarIH[i] = nuOsc->CalcMatterOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
      else
	{
	  t2kOscProbIH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMu, kPdgNuE, energy);
	  t2kOscProbBarIH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13+sinsq2theta13Error, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorS13IHPlus[i] = fabs(oscProb - t2kOscProbIH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13-sinsq2theta13Error, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorS13IHMinus[i] = fabs(oscProb - t2kOscProbIH[i]);

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23+sinsq2theta23PlusError, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorS23IHPlus[i] = fabs(oscProb - t2kOscProbIH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23-sinsq2theta23MinusError, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorS23IHMinus[i] = fabs(oscProb - t2kOscProbIH[i]);

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32-modDmsq32PlusError, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorDmsq32IHPlus[i] = fabs(oscProb - t2kOscProbIH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32+modDmsq32MinusError, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      t2kErrorDmsq32IHMinus[i] = fabs(oscProb - t2kOscProbIH[i]);

    }

  //---------------- Calculate P(numu->nue) and P(numubar->nuebar) for different values of deltaCP for NOvA -----------------------

  nuOsc->SetBaseline(810.0); 
  energy = 2.0;

  nuOsc->SetMassHierarchy(kNmhNormal);

  for(int i=0; i<npoints; i++)
    {
      deltaCP = i * 6.283185 / npoints;
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32, deltaCP);
      if(useMatterProb)
	{
	  novaNumuToNueNH[i] = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
	  novaNumubarToNuebarNH[i] = nuOsc->CalcMatterOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
      else
	{
	  novaNumuToNueNH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMu, kPdgNuE, energy);
	  novaNumubarToNuebarNH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
    }

  for(int i=0; i<4; i++)
    {
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      if(useMatterProb)
	{
	  novaOscProbNH[i] = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
	  novaOscProbBarNH[i] = nuOsc->CalcMatterOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
      else
	{
	  novaOscProbNH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMu, kPdgNuE, energy);
	  novaOscProbBarNH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13+sinsq2theta13Error, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorS13NHPlus[i] = fabs(oscProb - novaOscProbNH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13-sinsq2theta13Error, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorS13NHMinus[i] = fabs(oscProb - novaOscProbNH[i]);

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23+sinsq2theta23PlusError, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorS23NHPlus[i] = fabs(oscProb - novaOscProbNH[i]);;
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23-sinsq2theta23MinusError, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorS23NHMinus[i] = fabs(oscProb - novaOscProbNH[i]);

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32+modDmsq32PlusError, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorDmsq32NHPlus[i] = fabs(oscProb - novaOscProbNH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32-modDmsq32MinusError, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorDmsq32NHMinus[i] = fabs(oscProb - novaOscProbNH[i]);

    }

  nuOsc->SetMassHierarchy(kNmhInverted);

  for(int i=0; i<npoints; i++)
    {
      deltaCP = i * 6.283185 / npoints;
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32, deltaCP);
      if(useMatterProb)
	{ 
	  novaNumuToNueIH[i] = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
	  novaNumubarToNuebarIH[i] = nuOsc->CalcMatterOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
      else
	{
	  novaNumuToNueIH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMu, kPdgNuE, energy);
	  novaNumubarToNuebarIH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	} 
    }

  for(int i=0; i<4; i++)
    {
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      if(useMatterProb)
	{
	  novaOscProbIH[i] = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
	  novaOscProbBarIH[i] = nuOsc->CalcMatterOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}
      else
	{
	  novaOscProbIH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMu, kPdgNuE, energy);
	  novaOscProbBarIH[i] = nuOsc->CalcVacuumOscProb(kPdgNuMuBar, kPdgNuEBar, energy);
	}

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13+sinsq2theta13Error, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorS13IHPlus[i] = fabs(oscProb - novaOscProbIH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13-sinsq2theta13Error, sinsq2theta23, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorS13IHMinus[i] = fabs(oscProb - novaOscProbIH[i]);

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23+sinsq2theta23PlusError, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorS23IHPlus[i] = fabs(oscProb - novaOscProbIH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23-sinsq2theta23MinusError, dmsq12, modDmsq32, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorS23IHMinus[i] = fabs(oscProb - novaOscProbIH[i]);

      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32-modDmsq32PlusError, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorDmsq32IHPlus[i] = fabs(oscProb - novaOscProbIH[i]);
      nuOsc->SetOscParameters(sinsq2theta12, sinsq2theta13, sinsq2theta23, dmsq12, modDmsq32+modDmsq32MinusError, deltaCPpoints[i]);
      oscProb = nuOsc->CalcMatterOscProb(kPdgNuMu, kPdgNuE, energy);
      novaErrorDmsq32IHMinus[i] = fabs(oscProb - novaOscProbIH[i]);

    }

  //Plot graph on slide 3 (if useMatterProb is true) or graph on slide 5 (if useMatterProb is false)

  TGraph *T2KNH = new TGraph(npoints, t2kNumuToNueNH, t2kNumubarToNuebarNH);
  TGraph *T2KIH = new TGraph(npoints, t2kNumuToNueIH, t2kNumubarToNuebarIH);
  TGraph *NOvANH = new TGraph(npoints, novaNumuToNueNH, novaNumubarToNuebarNH);
  TGraph *NOvAIH = new TGraph(npoints, novaNumuToNueIH, novaNumubarToNuebarIH);

  TMarker *T2KNH00 = new TMarker(t2kOscProbNH[0], t2kOscProbBarNH[0], 25);
  TMarker *T2KNH05 = new TMarker(t2kOscProbNH[1], t2kOscProbBarNH[1], 26);
  TMarker *T2KNH10 = new TMarker(t2kOscProbNH[2], t2kOscProbBarNH[2], 27);
  TMarker *T2KNH15 = new TMarker(t2kOscProbNH[3], t2kOscProbBarNH[3], 30);
  TMarker *T2KIH00 = new TMarker(t2kOscProbIH[0], t2kOscProbBarIH[0], 25); 
  TMarker *T2KIH05 = new TMarker(t2kOscProbIH[1], t2kOscProbBarIH[1], 26);
  TMarker *T2KIH10 = new TMarker(t2kOscProbIH[2], t2kOscProbBarIH[2], 27);
  TMarker *T2KIH15 = new TMarker(t2kOscProbIH[3], t2kOscProbBarIH[3], 30);
  TMarker *NOvANH00 = new TMarker(novaOscProbNH[0], novaOscProbBarNH[0], 25);
  TMarker *NOvANH05 = new TMarker(novaOscProbNH[1], novaOscProbBarNH[1], 26);
  TMarker *NOvANH10 = new TMarker(novaOscProbNH[2], novaOscProbBarNH[2], 27);
  TMarker *NOvANH15 = new TMarker(novaOscProbNH[3], novaOscProbBarNH[3], 30);
  TMarker *NOvAIH00 = new TMarker(novaOscProbIH[0], novaOscProbBarIH[0], 25);
  TMarker *NOvAIH05 = new TMarker(novaOscProbIH[1], novaOscProbBarIH[1], 26);
  TMarker *NOvAIH10 = new TMarker(novaOscProbIH[2], novaOscProbBarIH[2], 27);
  TMarker *NOvAIH15 = new TMarker(novaOscProbIH[3], novaOscProbBarIH[3], 30);

  TCanvas *T2KNOvAComp = new TCanvas("T2KNOvAComp", "T2KNOvAComp", 800, 600);
  T2KNOvAComp->Divide(1);
  T2KNOvAComp->cd(1);
  FormatGraphLong(T2KNH, 2, 2, 1);
  T2KNH->GetXaxis()->SetTitle("P(#nu_{#mu} #rightarrow #nu_{e})");
  T2KNH->GetYaxis()->SetTitle("P(#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e})");
  T2KNH->SetLineColor(2);
  T2KNH->Draw("AL");
  FormatGraph(T2KIH, 4, 2, 1);
  T2KIH->Draw("L");
  FormatGraph(NOvANH, 8, 2, 1);
  NOvANH->Draw("L");
  FormatGraph(NOvAIH, 6, 2, 1);
  NOvAIH->Draw("L");
  T2KNH00->SetMarkerSize(2);
  T2KNH05->SetMarkerSize(2);
  T2KNH10->SetMarkerSize(3);
  T2KNH15->SetMarkerSize(2);
  T2KNH00->Draw();
  T2KNH05->Draw();
  T2KNH10->Draw();
  T2KNH15->Draw();
  T2KIH00->SetMarkerSize(2);
  T2KIH05->SetMarkerSize(2);
  T2KIH10->SetMarkerSize(3);
  T2KIH15->SetMarkerSize(2);
  T2KIH00->Draw();
  T2KIH05->Draw();
  T2KIH10->Draw();
  T2KIH15->Draw();
  NOvANH00->SetMarkerSize(2);
  NOvANH05->SetMarkerSize(2);
  NOvANH10->SetMarkerSize(3);
  NOvANH15->SetMarkerSize(2);
  NOvANH00->Draw();
  NOvANH05->Draw();
  NOvANH10->Draw();
  NOvANH15->Draw();
  NOvAIH00->SetMarkerSize(2);
  NOvAIH05->SetMarkerSize(2);
  NOvAIH10->SetMarkerSize(3);
  NOvAIH15->SetMarkerSize(2);
  NOvAIH00->Draw();
  NOvAIH05->Draw();
  NOvAIH10->Draw();
  NOvAIH15->Draw();

  TLegend *T2KNOvA = new TLegend(0.68, 0.48, 0.85, 0.87);
  T2KNOvA->AddEntry(T2KNH, "T2K NH", "L");
  T2KNOvA->AddEntry(T2KIH, "T2K IH", "L");
  T2KNOvA->AddEntry(NOvANH, "NOvA NH", "L");
  T2KNOvA->AddEntry(NOvAIH, "NOvA IH", "L");
  T2KNOvA->AddEntry(T2KNH00, "#delta_{CP} = 0", "P");
  T2KNOvA->AddEntry(T2KNH05, "#delta_{CP} = #pi/2", "P");
  T2KNOvA->AddEntry(T2KNH10, "#delta_{CP} = #pi", "P");
  T2KNOvA->AddEntry(T2KNH15, "#delta_{CP} = 3#pi/2", "P");
  T2KNOvA->SetBorderSize(0);
  T2KNOvA->Draw();

  //T2KNOvAComp->SaveAs("DeltaCPT2KNOvAComp.png");

  //Plot graph on slide 6 (if useMatterProb is true) or graph on slide 8 (if useMatterProb is false)

  TGraph *T2KNOvANH = new TGraph(npoints, t2kNumuToNueNH, novaNumuToNueNH);
  TGraph *T2KNOvAIH = new TGraph(npoints, t2kNumuToNueIH, novaNumuToNueIH);

  TMarker *T2KNOvANH00 = new TMarker(t2kOscProbNH[0], novaOscProbNH[0], 25);
  TMarker *T2KNOvANH05 = new TMarker(t2kOscProbNH[1], novaOscProbNH[1], 26);
  TMarker *T2KNOvANH10 = new TMarker(t2kOscProbNH[2], novaOscProbNH[2], 27);
  TMarker *T2KNOvANH15 = new TMarker(t2kOscProbNH[3], novaOscProbNH[3], 30);
  TMarker *T2KNOvAIH00 = new TMarker(t2kOscProbIH[0], novaOscProbIH[0], 25);
  TMarker *T2KNOvAIH05 = new TMarker(t2kOscProbIH[1], novaOscProbIH[1], 26);
  TMarker *T2KNOvAIH10 = new TMarker(t2kOscProbIH[2], novaOscProbIH[2], 27);
  TMarker *T2KNOvAIH15 = new TMarker(t2kOscProbIH[3], novaOscProbIH[3], 30);


  TCanvas *T2KNOvAComp2 = new TCanvas("T2KNOvAComp2", "T2KNOvAComp2", 800, 600);
  T2KNOvAComp2->Divide(1);
  T2KNOvAComp2->cd(1);
  FormatGraphLong(T2KNOvANH, 2, 2, 1);
  T2KNOvANH->GetXaxis()->SetTitle("T2K P(#nu_{#mu} #rightarrow #nu_{e})");
  T2KNOvANH->GetYaxis()->SetTitle("NOvA  P(#nu_{#mu} #rightarrow #nu_{e})");
  T2KNOvANH->Draw("AL");
  FormatGraph(T2KNOvAIH, 4, 2, 1);
  T2KNOvAIH->Draw("L");

  T2KNOvANH00->SetMarkerSize(2);
  T2KNOvANH05->SetMarkerSize(2);
  T2KNOvANH10->SetMarkerSize(3);
  T2KNOvANH15->SetMarkerSize(2);
  T2KNOvANH00->Draw();
  T2KNOvANH05->Draw();
  T2KNOvANH10->Draw();
  T2KNOvANH15->Draw();
  T2KNOvAIH00->SetMarkerSize(2);
  T2KNOvAIH05->SetMarkerSize(2);
  T2KNOvAIH10->SetMarkerSize(3);
  T2KNOvAIH15->SetMarkerSize(2);
  T2KNOvAIH00->Draw();
  T2KNOvAIH05->Draw();
  T2KNOvAIH10->Draw();
  T2KNOvAIH15->Draw();

  TLegend *T2KNOvA2 = new TLegend(0.22, 0.68, 0.57, 0.87);
  T2KNOvA2->AddEntry(T2KNOvANH, "Normal hierarchy", "L"); 
  T2KNOvA2->AddEntry(T2KNOvAIH, "Inverted hierarchy", "L");
  T2KNOvA2->SetBorderSize(0);
  T2KNOvA2->Draw();

  TLegend *T2KNOvA3 = new TLegend(0.66, 0.2, 0.87, 0.47);
  T2KNOvA3->AddEntry(T2KNOvANH00, "#delta_{CP} = 0", "P"); 
  T2KNOvA3->AddEntry(T2KNOvANH05, "#delta_{CP} = #pi/2", "P");
  T2KNOvA3->AddEntry(T2KNOvANH10, "#delta_{CP} = #pi", "P");
  T2KNOvA3->AddEntry(T2KNOvANH15, "#delta_{CP} = 3#pi/2", "P");
  T2KNOvA3->SetBorderSize(0);
  T2KNOvA3->Draw();
  
  //T2KNOvAComp2->SaveAs("HierarchyT2KNOvAComp.png");


  //------------------- Plot graphs with errors on slides 9 and 10 ------------------------------------------
  //---------------- (these are only plotted if useMatterProb is true)---------------------------------------

  if(useMatterProb)
    {   
    TGraphAsymmErrors *gS13ErrorNH = new TGraphAsymmErrors(4, t2kOscProbNH, novaOscProbNH, t2kErrorS13NHMinus,
							 t2kErrorS13NHPlus, novaErrorS13NHMinus, novaErrorS13NHPlus);
    TGraphAsymmErrors *gS23ErrorNH = new TGraphAsymmErrors(4, t2kOscProbNH, novaOscProbNH, t2kErrorS23NHMinus,
							 t2kErrorS23NHPlus, novaErrorS23NHMinus, novaErrorS23NHPlus);
    TGraphAsymmErrors *gDmsq32ErrorNH = new TGraphAsymmErrors(4, t2kOscProbNH, novaOscProbNH, t2kErrorDmsq32NHMinus,
							    t2kErrorDmsq32NHPlus, novaErrorDmsq32NHMinus, novaErrorDmsq32NHPlus);
    TGraphAsymmErrors *gS13ErrorIH = new TGraphAsymmErrors(4, t2kOscProbIH, novaOscProbIH, t2kErrorS13IHMinus,
							 t2kErrorS13IHPlus, novaErrorS13IHMinus, novaErrorS13IHPlus);
    TGraphAsymmErrors *gS23ErrorIH = new TGraphAsymmErrors(4, t2kOscProbIH, novaOscProbIH, t2kErrorS23IHMinus,
							 t2kErrorS23IHPlus, novaErrorS23IHMinus, novaErrorS23IHPlus);
    TGraphAsymmErrors *gDmsq32ErrorIH = new TGraphAsymmErrors(4, t2kOscProbIH, novaOscProbIH, t2kErrorDmsq32IHMinus,
							    t2kErrorDmsq32IHPlus, novaErrorDmsq32IHMinus, novaErrorDmsq32IHPlus);

    TCanvas *T2KNOvACompErrors = new TCanvas("T2KNOvACompErrors", "T2KNOvACompErrors", 800, 600);
    T2KNOvACompErrors->Divide(2, 2);
    FormatGraphLong(T2KNOvANH, 2, 2, 1);  
    T2KNOvANH->GetXaxis()->SetTitle("T2K P(#nu_{#mu} #rightarrow #nu_{e})");
    T2KNOvANH->GetYaxis()->SetTitle("NOvA  P(#nu_{#mu} #rightarrow #nu_{e})");
    T2KNOvANH->GetXaxis()->SetLimits(0.011, 0.089);
    T2KNOvANH->GetYaxis()->SetRangeUser(0.011, 0.089);
    FormatGraph(T2KNOvAIH, 4, 2, 1);  

    T2KNOvANH00->SetMarkerSize(2);
    T2KNOvANH05->SetMarkerSize(2);
    T2KNOvANH10->SetMarkerSize(3);
    T2KNOvANH15->SetMarkerSize(2);
    T2KNOvAIH00->SetMarkerSize(2);
    T2KNOvAIH05->SetMarkerSize(2);
    T2KNOvAIH10->SetMarkerSize(3);
    T2KNOvAIH15->SetMarkerSize(2);

    gS13ErrorNH->SetMarkerStyle(1);
    gS23ErrorNH->SetMarkerStyle(1);
    gDmsq32ErrorNH->SetMarkerStyle(1);
    gS13ErrorIH->SetMarkerStyle(1);
    gS23ErrorIH->SetMarkerStyle(1);
    gDmsq32ErrorIH->SetMarkerStyle(1);

    for(int i=0; i<3; i++)
      {
      T2KNOvACompErrors->cd(i+1);
      T2KNOvANH->Draw("AL");
      T2KNOvAIH->Draw("L");
      T2KNOvANH00->Draw();
      T2KNOvANH05->Draw();
      T2KNOvANH10->Draw();
      T2KNOvANH15->Draw();
      T2KNOvAIH00->Draw();
      T2KNOvAIH05->Draw();
      T2KNOvAIH10->Draw();
      T2KNOvAIH15->Draw();
      if(i == 0)
	{ 
	  gS13ErrorNH->Draw("P");
	  gS13ErrorIH->Draw("P");
	}
      else if(i == 1)
	{
	  gS23ErrorNH->Draw("P");
	  gS23ErrorIH->Draw("P");
	}
      else if(i == 2)
	{
	  gDmsq32ErrorNH->Draw("P");
	  gDmsq32ErrorIH->Draw("P");
	}

      TLegend *T2KNOvA4 = new TLegend(0.22, 0.68, 0.57, 0.87);
      T2KNOvA4->AddEntry(T2KNOvANH, "Normal hierarchy", "L");
      T2KNOvA4->AddEntry(T2KNOvAIH, "Inverted hierarchy", "L");
      T2KNOvA4->SetBorderSize(0);
      T2KNOvA4->Draw();

      TLegend *T2KNOvA5 = new TLegend(0.66, 0.2, 0.87, 0.47);
      T2KNOvA5->AddEntry(T2KNOvANH00, "#delta_{CP} = 0", "P");
      T2KNOvA5->AddEntry(T2KNOvANH05, "#delta_{CP} = #pi/2", "P");
      T2KNOvA5->AddEntry(T2KNOvANH10, "#delta_{CP} = #pi", "P");
      T2KNOvA5->AddEntry(T2KNOvANH15, "#delta_{CP} = 3#pi/2", "P");
      T2KNOvA5->SetBorderSize(0);
      T2KNOvA5->Draw();
      }

  //T2KNOvACompErrors->SaveAs("HierarchyT2KNOvACompErrors.png");

    }//end of if(useMatterProb)  

  //-------------------- Plot graph of comparison for antineutrinos on backup slide 19 ----------------------

  TGraph *T2KNOvABarNH = new TGraph(npoints, t2kNumubarToNuebarNH, novaNumubarToNuebarNH);
  TGraph *T2KNOvABarIH = new TGraph(npoints, t2kNumubarToNuebarIH, novaNumubarToNuebarIH);

  TMarker *T2KNOvABarNH00 = new TMarker(t2kOscProbBarNH[0], novaOscProbBarNH[0], 25);
  TMarker *T2KNOvABarNH05 = new TMarker(t2kOscProbBarNH[1], novaOscProbBarNH[1], 26);
  TMarker *T2KNOvABarNH10 = new TMarker(t2kOscProbBarNH[2], novaOscProbBarNH[2], 27);
  TMarker *T2KNOvABarNH15 = new TMarker(t2kOscProbBarNH[3], novaOscProbBarNH[3], 30);
  TMarker *T2KNOvABarIH00 = new TMarker(t2kOscProbBarIH[0], novaOscProbBarIH[0], 25);
  TMarker *T2KNOvABarIH05 = new TMarker(t2kOscProbBarIH[1], novaOscProbBarIH[1], 26);
  TMarker *T2KNOvABarIH10 = new TMarker(t2kOscProbBarIH[2], novaOscProbBarIH[2], 27);
  TMarker *T2KNOvABarIH15 = new TMarker(t2kOscProbBarIH[3], novaOscProbBarIH[3], 30);

  TCanvas *T2KNOvAComp3 = new TCanvas("HierarchyT2KNOvAComp3", "T2KNOvAComp3", 800, 600);
  T2KNOvAComp3->Divide(1);
  T2KNOvAComp3->cd(1);
  FormatGraphLong(T2KNOvABarNH, 2, 2, 1);
  T2KNOvABarNH->GetXaxis()->SetTitle("T2K P(#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e})");
  T2KNOvABarNH->GetYaxis()->SetTitle("NOvA  P(#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e})");
  T2KNOvABarNH->Draw("AL");
  FormatGraph(T2KNOvABarIH, 4, 2, 1);
  T2KNOvABarIH->Draw("L");

  T2KNOvABarNH00->SetMarkerSize(2);
  T2KNOvABarNH05->SetMarkerSize(2);
  T2KNOvABarNH10->SetMarkerSize(3);
  T2KNOvABarNH15->SetMarkerSize(2);
  T2KNOvABarNH00->Draw();
  T2KNOvABarNH05->Draw();
  T2KNOvABarNH10->Draw();
  T2KNOvABarNH15->Draw();
  T2KNOvABarIH00->SetMarkerSize(2);
  T2KNOvABarIH05->SetMarkerSize(2);
  T2KNOvABarIH10->SetMarkerSize(3);
  T2KNOvABarIH15->SetMarkerSize(2);
  T2KNOvABarIH00->Draw();
  T2KNOvABarIH05->Draw();
  T2KNOvABarIH10->Draw();
  T2KNOvABarIH15->Draw();

  TLegend *T2KNovA6 = new TLegend(0.22, 0.68, 0.57, 0.87);
  T2KNovA6->AddEntry(T2KNOvABarNH, "Normal hierarchy", "L");
  T2KNovA6->AddEntry(T2KNOvABarIH, "Inverted hierarchy", "L");
  T2KNovA6->SetBorderSize(0);
  T2KNovA6->Draw();

  TLegend *T2KNovA7 = new TLegend(0.67, 0.2, 0.87, 0.47);
  T2KNovA7->AddEntry(T2KNOvABarNH00, "#delta_{CP} = 0", "P");
  T2KNovA7->AddEntry(T2KNOvABarNH05, "#delta_{CP} = #pi/2", "P");
  T2KNovA7->AddEntry(T2KNOvABarNH10, "#delta_{CP} = #pi", "P");
  T2KNovA7->AddEntry(T2KNOvABarNH15, "#delta_{CP} = 3#pi/2", "P");
  T2KNovA7->SetBorderSize(0);
  T2KNovA7->Draw();

  //T2KNOvAComp3->SaveAs("HierarchyT2KNOvANubarComp.png");


}//end of plot_deltacp_hierarchy_comp_t2k_nova()


void FormatGraphLong(TGraph *graph, int colour, int width, int style)
{
  graph->SetTitle("");
  graph->SetLineColor(colour);
  graph->SetLineWidth(width);
  graph->SetLineStyle(style);
  graph->GetXaxis()->SetLimits(0.02, 0.08);
  graph->GetYaxis()->SetRangeUser(0.020001, 0.08);
  graph->GetXaxis()->SetTitleSize(0.06);
  graph->GetYaxis()->SetTitleSize(0.06);
  graph->GetXaxis()->SetLabelSize(0.06);
  graph->GetYaxis()->SetLabelSize(0.06);
  graph->GetXaxis()->SetTitleOffset(1.15);
  graph->GetYaxis()->SetTitleOffset(1.4);
  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->CenterTitle();
}

void FormatGraph(TGraph *graph, int colour, int width, int style)
{
  graph->SetTitle("");
  graph->SetLineColor(colour);
  graph->SetLineWidth(width);
  graph->SetLineStyle(style);
}
