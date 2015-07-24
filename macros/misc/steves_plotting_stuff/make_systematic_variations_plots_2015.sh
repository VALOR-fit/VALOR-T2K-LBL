
root2k -q -b ${T2KOSCANA}/macros/spectra/test_systematic_error_propagation_all_2015.C'('\
'0,false,"RHC1Rmu","3f_osc",'\
'0.304,0.0422,0.524,1.91,7.53E-5,2.51E-3,false,'\
'"analysis_configuration_2.5E20_prefit.xml","nuosc3f_12fix+nuisance","sinsqtheta","normal",'\
'-999,"f_banff_17_rhc",true,true)'

texfile='systematic_tables-sinsqtheta[12_0.304,13_0.042,23_0.524],dcp_1.910,dmsq_Normal[12_0.00007530,23_0.002510]eV2-RHC1Rmu'
if [ -f ${texfile}.tex  ] ; then
  pdflatex "${texfile}.tex"
  rm "${texfile}".{aux,log}
fi;

