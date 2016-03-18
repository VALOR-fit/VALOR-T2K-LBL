
double F1 = 0.5;
double F2 = 0.5;
double MA = 1.2257; // Axial Mass (BANFF post-fit)
double gA = 1.2694; // FA(Q2=0)
double Mp = 0.938;
double Eb = 0.027;
double Mlep = 0.000511;
double Mn = 0.939;

double CalcENu(double Plep, double costheta)
{
    double Elep = TMath::Power(Plep*Plep + Mlep*Mlep,0.5);
    
    return (Mp**2 - TMath::Power(Mn - Eb,2) - Mlep*Mlep + 2*(Mn - Eb)*Elep )/ (2*(Mn - Eb - Elep + Plep*costheta));
}

double CalcQ2(double Enu,double Elep,double Plep,double costheta)
{
    // Enu      --> neutrino reco energy
    // Elep     --> lepton reco energy
    // Plep     --> lepton reco momentum
    // costheta --> angle between neutrino and lepton
    
    // Q2 in the laboratory reference system (Giunti section 5.3.1)
    double Q2 = 2*Enu*(Elep-Plep*costheta) - Mlep*Mlep;
    
    return Q2;
}

double CalcXSecTerm(double* plep_1,double* costheta_1){
    
    // Term of CCQE dsigma/dQ^2 which sign is different for nu/nubar (arXiv:1305.7513v1, eq. 57-61)
    
    // Enu      --> neutrino reco energy
    // Elep     --> lepton reco energy
    // Plep     --> lepton reco momentum
    // costheta --> angle between neutrino and lepton
    // Mn       --> nucleon mass
    double Plep = plep_1[0];
    double costheta = plep_1[1];
    double Elep = TMath::Power(Plep*Plep + Mlep*Mlep,0.5);
    double Enu = CalcENu(Plep,costheta);

    double Q2 = CalcQ2(Enu,Elep,Plep,costheta);
    

    double FA = gA / TMath::Power( 1 + Q2/(MA*MA), 2. );
    double B = (Q2/(Mn*Mn))*FA*(F1+F2);
    
    double XSecTerm = ( B* (4*Mn*Enu - Q2 - Mlep*Mlep)) / (Mn*Mn);
    
    return XSecTerm;
}

double CalcFADerivative(double Q2, double costheta, double Q2Derivative)
{
    
    return (- 2* Q2Derivative * gA )/ ((MA*MA) * TMath::Power(1 + Q2/(MA*MA),3));
    
}

double CalcBDerivative(double Q2, double FA, double FADerivative, double Q2Derivative)
{
    
    return (Q2Derivative*FA*(F1+F2) - Q2 * FADerivative * (F1+F2) )/ (Mn*Mn) ;
}

double CalcQ2Derivative(double Enu,double Elep,double Plep,double costheta, double ENuDerivative)
{
    return 2* ENuDerivative*(Elep - Plep*costheta) + 2*Enu*(Plep / Elep - costheta);
}



double CalcENuDerivative(double Plep, double Elep, double costheta)
{
    return ((Mn - Eb)*(Plep/Elep) / (Mn - Eb -Elep +Plep*costheta))
    - ((Mn - Eb - (Plep/Elep) + costheta) *( 2*(Mn - Eb)* Elep + Mp*Mp - (Mn -Eb)**2 - Mlep*Mlep) )/ (2*TMath::Power(Mn - Eb -Elep +Plep*costheta,2));
}

double CalcXSecTermDerivative(double * par1, double* par2 )
{
    
    double Plep = par1[0];
    double costheta = par1[1];

    double Elep = TMath::Power(Mlep*Mlep + Plep*Plep,0.5);
    
    
    double Enu = CalcENu(Plep,costheta);
    
    double Q2 = CalcQ2(Enu,Elep,Plep,costheta);
    
    double FA = gA / TMath::Power( 1 + Q2/(MA*MA), 2. );
    double B = (Q2*FA*(F1+F2))/(Mn*Mn) ;
    
    
    double ENuDerivative = CalcENuDerivative(Plep,Elep,costheta);
    double Q2Derivative = CalcQ2Derivative(Enu,Elep,Plep,costheta,ENuDerivative);
    
    double FADerivative = CalcFADerivative(Q2,costheta,Q2Derivative );
    double BDerivative = CalcBDerivative(Q2,FA,FADerivative,Q2Derivative);
    
    double XSecDerivative = (((4*Mn*ENuDerivative - Q2Derivative) * B)  + ((4*Mn*Enu - Q2 - Mlep*Mlep) * BDerivative ) )/ (Mn*Mn);
//    
//cout << "Plep: " << Plep
//     << "\nCosT: " << costheta
//     << "\nElep: " << Elep
//     << "\nEnut: " << Enu
//     << "\nQsqd: " << Q2
//     << "\nFApr: " << FA
//     << "\nBpar: " << B
//     << "\ndEnu: " << ENuDerivative
//     << "\ndQsd: " << Q2Derivative
//     << "\ndFAd: " << FADerivative
//     << "\ndBpr: " << BDerivative
//     << "\ndXsc: " << XSecDerivative
//    << "\n\n" << endl;
    
    return XSecDerivative * (0.006*Plep + 0.026*TMath::Power(Plep,0.5));
}




XSec_term_and_error()
{
    TCanvas* c1 = new TCanvas("c","Error",800,500);

    TF2* xsecterm_error = new TF2("func",&CalcXSecTermDerivative,0,2, -1,1,2);
    TH2D* h = //new TH2D("h1","h1",30,0,2,20,-1,1);
    xsecterm_error->SetNpx(100);
    xsecterm_error->SetNpy(100);

    h = (TH2D*) xsecterm_error->CreateHistogram();
    h->Draw("colz");
    
    TCanvas* c2 = new TCanvas("c1","Term",800,500);
    TF2* xsecterm = new TF2("xsecterm",&CalcXSecTerm,0,2, -1,1,2);
    xsecterm->SetNpx(100);
    xsecterm->SetNpy(100);

    TH2D* h2 = //new TH2D("h1","h1",30,0,2,20,-1,1);
    h2 = (TH2D*) xsecterm->CreateHistogram();
    h2->Draw("colz");
    
}
