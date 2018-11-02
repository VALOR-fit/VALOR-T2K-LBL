# T2K-LBL
VALOR T2K(/HK) LBL oscillation analysis

# T2K LBL analysis releases and associated framework releases

T2K-LBL tag | FMWK tag | Summary
-----------------------------------------------------------------------------------------------------------------------------
2010_v1.0.0 |   n/a    | This is the code which was used for the 2010a 0.323E+20 POT muon-disappearance analysis.
            |          | Details may be found in T2K-TN-36 (v4).
            |          |
2010_v2.0.0 |   n/a    | This is the code which was used for the 2010b 1.431E+20 POT muon-disappearance analysis.
            |          | Details may be found in T2K-TN-64 (v3).
            |          |
2010_v2.1.0 |   n/a    | Similar to 2010_v2.0.0 but with the following change requested during the analysis review:
            |          | > The nd280 efficiency error was increased from +/-5% to +5.6% -5.2% to match the one used 
            |          |   in the nue analysis.
            |          | > The statistical error of the ND280 CC inclusive measurement (~2.7%) was also propagated to SK.
            |          |
2010_v3.0.0 |   n/a    | This code was used for the 2010b 1.431E+20 POT 3-flavour muon-neutrino disappearance analysis 
            |          |(details in T2K-TN-88) and the 2010b 1.431E+20 POT joint 3-flavour analysis (details in T2K-TN-131)            
            |          |
2012_v1.0.0 |   n/a    | This is the code used for the official 2012 3+E_20 POT 3-flavour muon-neutrino disappearance analysis.
            |          | A *major* re-write of the analysis wrt to the 2010 releases since the analysis inputs and
            |          | the analysis procedure (fitting for systematics) were very different.
            |          | Details, too numerous to mention here, may be found in T2K-TN-141 (v4).
            |          | This is the first freeze of the fully validated 2012 code: This branch also contains experimental and 
            |          | deprecated features (experimental FSI and BANFF binnings, systematics not used in the final analysis,
            |          | deprecated macros etc) that will be removed from a near-future "streamlined" revision.
            |          |
2012_v1.1.0 |   n/a    | Revision with small bugs fixed while writing the 2012 analysis PRL and resolving final differences with KIT
            |          |
2012_v1.2.0 |   n/a    | Revision with all changes needed during the PRL review process. Single most important change was moving from
            |          | sin^2(2theta23) fit to sin^2(theta23) fit.
            |          |
2012_v2.0.0 |   n/a    | Implemented the 2012 3-flavour electron-neutrino appearance analysis and the 2012 joint 3-flavour oscillation 
            |          | analysis. The source tree was cleaned-up (Experimental FSI and BANFF systematic binnings, systematics not 
            |          | used in the analysis were removed.  Deprecated data files, macros and scripts were removed.  
            |          | Code experimenting with the use of principal  components analysis in the study of the effect of systematics 
            |          | in the oscillation fit was deleted.) The main Minuit and the experimental Minuit2 fit implementations were 
            |          | "unified" so that one can switch between the two  transparently. A GSL fitter was implementing, adding the 
            |          | following fit methods in our analysis: Simulated annealing, Levemberg-Marquardt, Fletcher-Reeves 
            |          | conjugate gradient, Polak-Ribiere conjugate gradient, Vector Broyden-Fletcher-Goldfarb-Shanno (BFGS), 
            |          | Steepest Descent. Significant fit CPU-efficiency improvement as the jont analysis considers ~90 systematics.
            |          |
2013_v1.0.0 |   n/a    | This is the code used for the official 2013 Run 1-4a 3-flavour muon-neutrino disappearance analysis.
            |          | Details may be found in T2K-TN-155v1.
            |          |
2013_v2.0.0 |   n/a    | This is the code used for the official 2013 Run 1-4 3-flavour muon-neutrino disappearance analysis.
            |          | Details may be found in T2K-TN-183v1.
            |          |
2013_v3.0.0 |   n/a    | This is the code used for the early versions of the 2013 Run 1-4 joint analysis.
            |          | Details may be found in T2K-TN-175
            |          | 
2013_v3.1.0 |   n/a    | This is the code used for the final versions of the 2013 Run 1-4 joint analysis.
            |          | Details may be found in T2K-TN-175
            |          |
2014_v1.0.0 |  v1.0.0  | First version of the code within the new VALOR framework. Not used for an analysis.
            |          |
2014_v1.1.0 |  v1.0.2  | Freeze before moving to the 2014/15 cross-section parameterisation. Not used for analysis.
            |          |
2015_v1.0.0 |  v1.1.0  | This is the code used for the 2015 Run 5+6 3-flavour muon-antineutrino disappearance analysis.
            |          | This version was used for the first version of the tech note (2.5E20 POT, data not yet fit)
            |          | Details may be found in T2K-TN-243v0.
            |          |
2015_v1.1.0 |  v1.1.0  | This is the code used for the 2015 Run 5+6 3-flavour muon-antineutrino disappearance analysis.
            |          | This version was used for the final version of the tech note (2.315E20 POT, including data fit)
            |          | Details may be found in T2K-TN-243v3.
            |          |
2015_v2.0.0 |  v1.1.0  | This is the code used for the 2015 Run 5+6 3-flavour electron antineutrino appearance analysis.
            |          | This version was used for the first version of the tech note (2.5E20 POT, data not yet fit)
            |          | Details may be found in T2K-TN-252v1.
            |          |
2016_v1.0.0 |  v1.2.0  | This is the code used for the Neutrino 2016 analyses
            |          |
            |          |
2017_v0.0.1 |  v1.4.2  | This is the code for the 2017 joint analysis (not yet frozen)
            |          |
            |          |
