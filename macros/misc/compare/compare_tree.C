//
// Compare 2 TTree's for compatibility
// Currently set up for comparing the 2 SKMC trees we use in the analysis
//  1 full tree with all variables used for reweighting (to generate splines, FSI/SI matrix)
//  2 summary tree with just Etrue, Ereco, NEUT mode
//

#include <stdlib.h>
#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>

#define fASmallNum 0.000001

bool Compare(int a, int b){
  bool arethesame = TMath::Abs(a-b) < fASmallNum;
  if(!arethesame){ cout << "Compared: "<< a <<" to " << b << " --> diff = "<< a-b << endl; }
  return arethesame;
}

bool Compare(unsigned a, unsigned b){
  bool arethesame = TMath::Abs((int)a-(int)b) < fASmallNum;
  if(!arethesame){ cout << "Compared: "<< a <<" to " << b << " --> diff = "<< a-b << endl; }
  return arethesame;
}

bool Compare(double a, double b){
  bool arethesame = TMath::Abs(a-b) < fASmallNum;
  if(!arethesame){ cout << "Compared: "<< a <<" to " << b << " --> diff = "<< a-b << endl; }
  return arethesame;
}

void compare_tree(
  const char * filename1 = "mtuple_sk11c_flx11bv31_mu10acuts_numu.root",               /* name of 1st file */
  const char * filename2 = "$T2KOSCANA/inputs/official/2012/superk/mc/neut/11c/smallnt_numu_sk11c_flx11bv31fine_r1-3cfull_10acut.root",               /* name of 2nd file */
  string treename1 = "mtuple",  /* name of 1st tree */
  string treename2 = "tnumu"                /* name of 2nd tree */
)

{

  //get trees
  TFile * f1 = new TFile(filename1,"read");
  TFile * f2 = new TFile(filename2,"read");

  if(treename2 == "")
    treename2 = treename1;

  TTree * t1 = (TTree*) f1->Get(treename1.c_str());
  TTree * t2 = (TTree*) f2->Get(treename2.c_str());

  if(t1->GetEntries() != t2->GetEntries()) {
    cout << t1->GetEntries() << "!=" << t2->GetEntries() << " number of entries don't match!" << endl;
    exit(1);
  }
  cout << "Number of entries is " <<  t1->GetEntries() << endl;

  int failed=0;

  //set tree branches
  int Mode1=0, Mode2=0;
  float Enu1[50], Enu2=0, EnuReco2=0;
  double EnuReco1=0;

  t1->SetBranchAddress("mode",&Mode1);
  t1->SetBranchAddress("pnu",Enu1);
  t1->SetBranchAddress("erec",&EnuReco1);
  
  t2->SetBranchAddress("mode",&Mode2);
  t2->SetBranchAddress("enu",&Enu2);
  t2->SetBranchAddress("erec",&EnuReco2);

  //loop over tree
  for(int i=0; i<t1->GetEntries(); i++) {
    t1->GetEntry(i);
    t2->GetEntry(i);
    
    if (!Compare(Mode1,Mode2)) {failed++; cout << "failed check for mode " << Mode1 << "," << Mode2 << endl;}
    if (!Compare(Enu1[0],Enu2)) {failed++; cout << "failed check for enu " << Enu1 << "," << Enu2 << endl;}
    if (!Compare((float)EnuReco1,EnuReco2)) {failed++; cout << "failed check for enuReco " << EnuReco1 << "," << EnuReco2 << endl;}
  }//loop over tree entries

  cout << "Total failures is: " << failed << endl;
}
