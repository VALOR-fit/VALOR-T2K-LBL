//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include <cassert>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <TSystem.h>

#include "Enumerations.h"
#include "Conventions.h"
#include "Utilities.h"

#include "UtlFitParamListFactory.h"
#include "FitCore/FitParamList.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace valor;
using namespace valor::conv;

FitParamList * util::NamedParameterList(string name)
{

  //
  // 3-flavour oscillation fit with/without systematics included 
  //

  if(name == "nuosc3f") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }

  //
  // 3-flavour oscillation fit with various oscillation param subsets either fixed or
  // treated as nuisance paramers.
  //

  if(name == "nuosc3f_12constrain") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f_12constrain+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_12fix") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f_12fix+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }

  if(name == "nuosc3f_numudisapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  
  if(name == "nuosc3f_numudisapp_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  
  if(name == "nuosc3f_numudisapp_nuis3fixdcp+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance (kPhypMixing13, true);
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq, true);
    fp -> FixParameter (kPhypCPV);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }

  if(name == "nuosc3f_nueapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f_nueapp_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix5") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix3") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(4);//use future osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(3);//use current osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix3+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(5);//use future osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(6);//use current osc errors
    return fp;
  }
  
  cerr<<"FitParamListFactory/NamedParamList >> Unable to make: "<<name<<endl;
  exit(-1); 
}

//________________________________________________________________________________________
void util::AppendAllSystematicsTo(FitParamList * fparm, int start)
{
// Add nuisance params
//

#ifdef __NUMU_ONLY__
  util::AppendAllCorrelatedBANFFSystematicsDisap  (fparm, start);
  util::AppendAllSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfParamsInBANFFfit);
  util::AppendAllUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
#elsif __NUE_ONLY__
  util::AppendAllCorrelatedBANFFSystematicsApp  (fparm, start);
  util::AppendAllSuperkDetectorFSISystematicsApp(fparm, start+kNumOfParamsInBANFFfit);
  util::AppendAllUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
#else
  util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
  util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfParamsInBANFFfit);
  util::AppendAllUncorrelatedXSecSystematics   (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
#endif

}
//________________________________________________________________________________________

//________________________________________________________________________________________
void util::AppendAllUncorrelatedXSecSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecFermiMom));
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecSpecFunc));
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecWshape));
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecShapeCCoth));
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecPiLessDelta));
#ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecBindE));
#endif
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormCCcoh));
#ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNC1piC));
#endif
#ifndef __NUMUONLY__
  fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP*0.3,  false, true, AsString(kNuisXSecNormNCcoh));
#endif
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCoth));
#ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormAntiNuToNu));
#endif
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNuEToNuMu));


}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void util::AppendUncorrelatedXSecSystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list
  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecFermiMom) );
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecSpecFunc) );
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecWshape) );
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecShapeCCoth) );
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecPiLessDelta) );
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecBindE) );
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormCCcoh) );
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNC1piC) );
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCoth) );
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNuEToNuMu) );
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormAntiNuToNu) );
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendUncorrelatedXSecSystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list
  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecFermiMom) );
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecSpecFunc) );
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecWshape) );
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecShapeCCoth) );
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecPiLessDelta) );
  //  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecBindE) );
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormCCcoh) );
  fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCcoh) );
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCoth) );
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNuEToNuMu) );
  //  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormAntiNuToNu) );
}
#endif
//________________________________________________________________________________________
void util::AppendAllSuperkDetectorFSISystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list
  // SK detector energy scale
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale) );

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get FSI covariance matrix
#if FSI_BINS == 1
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/covariance_matrix");
#elif FSI_BINS == 2
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/disap/covariance_matrix");
#elif FSI_BINS == 3
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/app/covariance_matrix");
#endif
  if(!SKDetFSICovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No SKDET+FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDetFSICov = *SKDetFSICovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // FSI
  for(int i = 0; i < kNumOfSKDetFSIParams; i++) {
#if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
#else
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
#endif
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void util::AppendSuperkDetectorFSISystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list
  // SK detector energy scale
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get FSI covariance matrix
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/disap/covariance_matrix");
  if(!SKDetFSICovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No SKDET+FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDetFSICov = *SKDetFSICovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // FSI
  for(int i = 0; i < kNumOfSKDetFSIParams1Rmu; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendSuperkDetectorFSISystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list
  // SK detector energy scale
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get FSI covariance matrix
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/app/covariance_matrix");
  if(!SKDetFSICovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No SKDET+FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDetFSICov = *SKDetFSICovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // FSI
  for(int i = 0; i < kNumOfSKDetFSIParams1Re; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
}
#endif
//________________________________________________________________________________________
void util::AppendAllCorrelatedBANFFSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get BANFF fit central values
#if BANFF_BINS == 1
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/best_fit_values");
#elif BANFF_BINS == 2
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/disap/best_fit_values");
#elif BANFF_BINS == 3
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/app/best_fit_values");
#endif
  if(!BANFFBestFitValues) {
    cerr << "FitParamList/AppendSystematics >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
#if BANFF_BINS == 1
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/covariance_matrix");
#elif BANFF_BINS == 2
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/disap/covariance_matrix");
#elif BANFF_BINS == 3
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/app/covariance_matrix");
#endif
  if(!BANFFCovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & BANFFCov = *BANFFCovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // BANFF 
  for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
    PhysParamId_t p = kListOfBANFFParams[i];
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }

}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void util::AppendCorrelatedBANFFSystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get BANFF fit central values
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/disap/best_fit_values");
  if(!BANFFBestFitValues) {
    cerr << "FitParamList/AppendSystematics >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/disap/covariance_matrix");
  if(!BANFFCovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & BANFFCov = *BANFFCovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // BANFF 
  for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
    PhysParamId_t p = kListOfBANFFParams[i];
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendCorrelatedBANFFSystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get BANFF fit central values
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/app/best_fit_values");
  if(!BANFFBestFitValues) {
    cerr << "FitParamList/AppendSystematics >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/app/covariance_matrix");
  if(!BANFFCovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & BANFFCov = *BANFFCovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // BANFF 
  for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
    PhysParamId_t p = kListOfBANFFParams[i];
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
}
#endif
//________________________________________________________________________________________
