#include <iostream>
#include <iomanip>

using namespace std;

void read_fitparamlist_obj(
   const char * filename = "./fp.root", const char * objname = "fit_param_list")
{
  using namespace valor;

  // open ROOT file
  TFile f(filename,"read");

  // get a param list
  FitParamList * fparm = (FitParamList *)f.Get(objname);
  assert(fparm);

  // print stuff
  cout << "Number of params: " << fparm->GetNumberParameters() << endl;
  cout << "Number of free params: " << fparm->GetNumberFreeParameters() << endl;

  cout << *fparm;

  // close
  f.Close();
}
