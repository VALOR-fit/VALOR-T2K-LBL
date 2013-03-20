/*
Takes the output of the star plot (in tree form), finds the maximum deviation from nominal (for a given parameter,
 as a fraction of the stat error, from a variety of points)

Tom
 */

void bestfit_point_shift_table(
 string       data_file = "star_plot_line-uncor_test2.data", /* merged outputs of test_bestfit_point_shift_due_to_syst_twk.C calls */
 bool         correlated = false,
 double       true_mixing23 = 0.4,
 double       true_dm23sq = 2.4E-3,
 double       true_mixing13 = 0.025,
 double       true_deltacp = 0,
 double       stat_av_mixing23 = 0.0605587,
 double       stat_av_dm23sq = 0.000181143,
 double       stat_av_mixing13 = 0.013045,
 double       stat_av_deltacp = 1.90991
)
{
  using namespace std;
  using namespace t2kana;


  ifstream fit_data(data_file.c_str());
  if (!fit_data.good() ) {
      cerr << "** Can't open data file: " << data_file << endl;
      return;
  }

  int nprecisionnev = 4;
  stringstream streamer;  

  streamer << "\\begin{table}"   << endl;
  streamer << "\\scriptsize" << endl;
  streamer << "\\begin{tabular} {|p{2cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}|}" << endl;
  streamer << "\\hline" << endl;
  streamer << " {Systematic}  & \\multicolumn {4}{c||} {+ 1 ${\\sigma}$ } & \\multicolumn {4}{c|} {- 1 ${\\sigma}$ } \\\\ " << endl; 
  streamer << "\\hline \\hline" << endl;
  streamer << "  &  {${sin^{2}(\\theta_{23})}$} & {$|{\\Delta m^{2}}|$} & {${sin^{2}(\\theta_{13})}$} & {${\\delta_{CP}}$}  & {${sin^{2}(\\theta_{23})}$} &  {$|{\\Delta m^{2}}|$} & {${sin^{2}(\\theta_{13})}$}  &  {${\\delta_{CP}}$}  \\\\ " << endl; 
  streamer << "\\hline" << endl;
  streamer << "  & {Deviation (\\%)} &  {Deviation (\\%)} & {Deviation (\\%)}& {Deviation (\\%)} & {Deviation (\\%)} & {Deviation (\\%)} & {Deviation (\\%)}  & {Deviation (\\%)}  \\\\ " << endl; 
  streamer << "\\hline" << endl;

  int isyst=0;
  while(1) {

      // skip header lines staring with #
      if(fit_data.peek() == '#') {
         fit_data.ignore(1000, '\n');

       } else {
	const char * systematic = "";
	double       bestfit_mixing13_p1sig  = 0.;
	double       bestfit_mixing23_p1sig  = 0.;
	double       bestfit_deltacp_p1sig   = 0.;
	double       bestfit_dm23sq_p1sig    = 0.;
	double       bestfit_mixing13_m1sig  = 0.;
	double       bestfit_mixing23_m1sig  = 0.;
	double       bestfit_deltacp_m1sig   = 0.;
	double       bestfit_dm23sq_m1sig    = 0.;
	
	fit_data >> systematic 
                  >> bestfit_mixing13_p1sig
                  >> bestfit_mixing23_p1sig
                  >> bestfit_deltacp_p1sig
                  >> bestfit_dm23sq_p1sig
                  >> bestfit_mixing13_m1sig
                  >> bestfit_mixing23_m1sig
                  >> bestfit_deltacp_m1sig
		 >> bestfit_dm23sq_m1sig;

	if(fit_data.eof()) break;

	double mixing13_dev_p1sig = (TMath::Abs(bestfit_mixing13_p1sig - true_mixing13)*100.)/stat_av_mixing13;
	double mixing23_dev_p1sig = (TMath::Abs(bestfit_mixing23_p1sig - true_mixing23)*100.)/stat_av_mixing23;
	double dm23sq_dev_p1sig   = (TMath::Abs(bestfit_dm23sq_p1sig - true_dm23sq)*100.)/stat_av_dm23sq;
	double deltacp_dev_p1sig  = (TMath::Abs(bestfit_deltacp_p1sig - true_deltacp)*100.)/stat_av_deltacp;

	double mixing13_dev_m1sig = (TMath::Abs(bestfit_mixing13_m1sig - true_mixing13)*100.)/stat_av_mixing13;
	double mixing23_dev_m1sig = (TMath::Abs(bestfit_mixing23_m1sig - true_mixing23)*100.)/stat_av_mixing23;
	double dm23sq_dev_m1sig   = (TMath::Abs(bestfit_dm23sq_m1sig - true_dm23sq)*100.)/stat_av_dm23sq;
	double deltacp_dev_m1sig  = (TMath::Abs(bestfit_deltacp_m1sig - true_deltacp)*100.)/stat_av_deltacp;


	streamer << " {" << util::AsStringPaperTex(util::PhysParamIdFromString(systematic)) << "}  & " ;
	if(mixing23_dev_p1sig > 25) streamer << "\\color{red}";
	streamer << "{ " << setprecision(nprecisionnev) << mixing23_dev_p1sig << "} & " ;
	if(dm23sq_dev_p1sig > 25) streamer << "\\color{red}";
	streamer << "{ " << setprecision(nprecisionnev) << dm23sq_dev_p1sig << "} & " ;
	if(mixing13_dev_p1sig > 25) streamer << "\\color{red}";
	streamer << "{ " << setprecision(nprecisionnev) << mixing13_dev_p1sig << "} & " ;
	if(deltacp_dev_p1sig > 25) streamer << "\\color{red}";
	streamer << "{ " << setprecision(nprecisionnev) << deltacp_dev_p1sig << "} & " ;

	if(mixing23_dev_m1sig > 25) streamer << "\\color{red}";
	streamer << "{ " << setprecision(nprecisionnev) << mixing23_dev_m1sig << "} & " ;
	if(dm23sq_dev_m1sig > 25) streamer << "\\color{red}";
	streamer << "{ " << setprecision(nprecisionnev) << dm23sq_dev_m1sig <<  "} & " ;
	if(mixing13_dev_m1sig > 25) streamer << "\\color{red}";
	streamer << "{ " << setprecision(nprecisionnev) << mixing13_dev_m1sig << "} & " ;
	if(deltacp_dev_m1sig > 25) streamer << "\\color{red}";
	streamer << "{ " << setprecision(nprecisionnev) << deltacp_dev_m1sig << "} \\\\" << endl;

	streamer << "\\hline" << endl;

	if (isyst == 44){ //divide into two tables: first one BANFF + Xsec, second eff + FSI
	  streamer << "\\end{tabular}" << endl;
	  streamer << "\\caption{Effect of BANFF and cross section systematics in the best fit point for true values: ${sin^{2}(\\theta_{13})}$ = " << true_mixing13 << " , ${sin^{2}(\\theta_{23})}$ =  " << true_mixing23 << " , $\\delta_{CP}$ = "<< true_deltacp << "  and {$|{\\Delta m^{2}}|$} = " << true_dm23sq << "} " << endl;
	  streamer << "\\end{table}" << endl;
	  streamer << "\\begin{table}"   << endl;
	  streamer << "\\scriptsize" << endl;
	  // new table
	  streamer << "\\begin{tabular} {|p{2cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}||p{1cm}|}" << endl;
	  streamer << "\\hline" << endl;
	  streamer << " {Systematic}  & \\multicolumn {4}{c||} {+ 1 ${\\sigma}$ } & \\multicolumn {4}{c|} {- 1 ${\\sigma}$ } \\\\ " << endl; 
	  streamer << "\\hline \\hline" << endl;
	  streamer << "  &  {${sin^{2}(\\theta_{23})}$} & {$|{\\Delta m^{2}}|$} & {${sin^{2}(\\theta_{13})}$} & {${\\delta_{CP}}$}  & {${sin^{2}(\\theta_{23})}$} &  {$|{\\Delta m^{2}}|$} & {${sin^{2}(\\theta_{13})}$}  &  {${\\delta_{CP}}$}  \\\\ " << endl; 
	  streamer << "\\hline" << endl;
	  streamer << "  & {Deviation (\\%)} &  {Deviation (\\%)} & {Deviation (\\%)}& {Deviation (\\%)} & {Deviation (\\%)} & {Deviation (\\%)} & {Deviation (\\%)}  & {Deviation (\\%)}  \\\\ " << endl; 
	  streamer << "\\hline" << endl;
	}
	isyst++;
      }//not a comment line
  }//while(1)

  streamer << "\\end{tabular}" << endl;
  streamer << "\\caption{Effect of SK efficiencies and FSI+SI systematics in the best fit point for true values: ${sin^{2}(\\theta_{13})}$ = " << true_mixing13 << " , ${sin^{2}(\\theta_{23})}$ =  " << true_mixing23 << " , $\\delta_{CP}$ = "<< true_deltacp << "  and {$|{\\Delta m^{2}}|$} = " << true_dm23sq << "} " << endl;
  streamer << "\\end{table}" << endl;
    

 //  cout << streamer.str() << endl;
  long size;
  // get size of file                                                                                                                 
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);
  
  ostringstream output_file;
  string cor = "";
  if(correlated)
    cor = "_cor_";
  else cor = "_uncor_";
  output_file << "Effect_systematics_best_fit_point" << cor << true_mixing23 << "_" << true_dm23sq << "_" << true_mixing13 << "_" << true_deltacp << ".tex";
  ofstream outfile (output_file.str().c_str(),ios_base::app);

  // write document header to outfile                  
  string f;
  f += "\\documentclass[10pt,a4paper]{article}\n";
  f += ("\\usepackage{amsmath}\n");
  f += ("\\usepackage{graphicx}\n");
  f += ("\\usepackage{dcolumn}\n");
  f += ("\\usepackage[cm]{fullpage}\n");
  f += ("\\usepackage{color}\n");
  f += ("\\begin{document}\n");
  outfile.write (f.c_str(),f.size());
  //outfile.close();
  
  // write to outfile                                  
  outfile.write (streamer.str().c_str(),size);
  
  //after loop over modes 
  string f;
  f += ("\\end{document}\n");
   outfile.write (f.c_str(),f.size());
  outfile.close();



}
