#include <iostream>
#include <iomanip>

using namespace std;

void read_physparamlist_obj(
   const char * filename = "./physlist.root", const char * objname = "phys_param_list")
{
  using namespace valor;

  // open ROOT file
  TFile f(filename,"read");

  // get the param list
  PhysParamList * plptr = (PhysParamList *)f.Get(objname);
  assert(plptr);

  PhysParamList & pl = *plptr;
  
  // prnt stuff
  cout << pl;

  cout << "There were " << pl.NParams() << " parameters in the list" << endl;
  cout << "Parameter " << util::AsString(conv::kPhypMixing23) << " value was " << pl(conv::kPhypMixing23) << endl;

  // close
  f.Close();
}
