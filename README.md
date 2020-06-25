# VALOR T2K-LBL

The VALOR-T2K-LBL repository includes T2K-specific codes and configurations for the VALOR T2K oscillation fit and related analysis workflows. This analysis is built on top of the common and experiment-agnostic [VALOR *Software Development Kit (SDK)*](https://github.com/VALOR-fit/VALOR-SDK/) which includes the bulk of the VALOR analysis code.

**Copyright (c) 2010-2020**, The VALOR Fitting Group. All Rights Reserved.

**Authors (VALOR/T2K):** 
* Costas Andreopoulos (*University of Liverpool & STFC/RAL*)
* Maria Antonova (*IFIC Valencia*)
* Chistopher Barry (*University of Liverpool*)
* Francis Bench (*University of Liverpool*)
* Andrew Chappell (*University of Warwick*)
* Τhomas Dealtry, (*Lancaster University; Oxford & STFC/RAL*)
* Steve Dennis (*University of Liverpool; University of Warwick & STFC/RAL*)
* Lorena Escudero (*IFIC Valencia*)
* Pablo Fernandez Menendez (*IFIC Valencia*)
* Jaiden Parlone (*IFIC Valencia*)
* Davide Sgalaberna (*CERN; Geneva; ETHZ*)
* Raj Shah (*University of Oxford & STFC/RAL*)

For a full list of VALOR SDK authors, see [here](https://valor.pp.rl.ac.uk/).


## T2K LBL oscillation analysis releases 

Since 2010, the VALOR group led numerous flagship T2K oscillation analyses, produced world-leading results and contributed to most published T2K oscillation papers to date. The version of the code that can be used for reproducing each analysis, can be found below" 

<pre>

VALOR-T2K-LBL| VALOR-SDK      | Summary
--------------------------------------------------------------------------------------------------------------------
2010_v1.0.0  |   n/a    | Code used for the 2010 Run 1 2-flavour muon-disappearance analysis.
             |          | Details may be found in T2K-TN-36 (v4).
             |          |
2010_v2.0.0  |   n/a    | Code used for the 2010 Run 1-2 2-flavour muon-disappearance analysis.
             |          | Details may be found in T2K-TN-64 (v1).
             |          |
2010_v2.1.0  |   n/a    | Similar to 2010_v2.0.0 but including changes requested during the review.
             |          | Details may be found in the final revision of T2K-TN-64 (v3).
             |          |
2010_v3.0.0  |   n/a    | Code used for the 2010 Run 1-2 3-flavour muon-disappearance analysis. 
             |          | Details may be found in T2K-TN-87.
             |          |
2012_v1.0.0  |   n/a    | Code used for the official 2012 Run 1-3 3-flavour muon-neutrino disappearance analysis.
             |          | Details may be found in T2K-TN-141 (v4).
             |          | A major re-write of the analysis code.
             |          | This is the first freeze of the fully validated 2012 code which contains experimental  
             |          | and deprecated features that will be removed from a near-future revision.
             |          |
2012_v1.1.0  |   n/a    | Revision with small bugs fixed while writing the 2012 analysis PRL 
             |          | and resolving final differences with KIT
             |          |
2012_v1.2.0  |   n/a    | Revision with all changes needed during the PRL review process. 
             |          | Single most important change was moving from sin^2(2theta23) fit to sin^2(theta23) fit.
             |          |
2012_v2.0.0  |   n/a    | First implementation of the 2012 Run 1-3 3-flavour electron-neutrino appearance analysis 
             |          | and of a joint 3-flavour oscillation analyses within the VALOR framework!
             |          | The source tree was cleaned-up and ignificant fit CPU-efficiency improvements were made.
             |          | Details may be found in T2K-TN-154.
             |          |
2013_v1.0.0  |   n/a    | Code used for the official 2013 Run 1-4a 3-flavour muon-neutrino disappearance analysis.
             |          | Details may be found in T2K-TN-155v1.
             |          |
2013_v2.0.0  |   n/a    | Code used for the official 2013 Run 1-4 3-flavour muon-neutrino disappearance analysis.
             |          | Details may be found in T2K-TN-183v1.
             |          |
2013_v3.0.0  |   n/a    | Code used for the early versions of the 2013 Run 1-4 joint analysis.
             |          | Details may be found in T2K-TN-175
             |          | 
2013_v3.1.0  |   n/a    | Code used for the final versions of the 2013 Run 1-4 joint analysis.
             |          | Details may be found in T2K-TN-175
             |          |
2014_v1.0.0  |  v1.0.0  | An early version of the T2K LBL analysis, factoring out an experiment-agnostic VALOR SDK.
             |          | Not used for an analysis.
             |          |
2014_v1.1.0  |  v1.0.2  | Freeze before moving to the 2014/15 cross-section parameterisation. 
             |          | Not used for analysis.
             |          |
2015_v1.0.0  |  v1.1.0  | Code used for the 2015 Run 5-6 3-flavour muon-antineutrino disappearance analysis.
             |          | This version was used for the first version of the tech note. 
             |          | (2.5E20 POT, data not yet fit)
             |          | Details may be found in T2K-TN-243v0.
             |          |
2015_v1.1.0  |  v1.1.0  | Code used for the 2015 Run 5-6 3-flavour muon-antineutrino disappearance analysis.
             |          | This version was used for the final version of the tech note. 
             |          | (2.315E20 POT, including data fit)
             |          | Details may be found in T2K-TN-243v3.
             |          |
2015_v2.0.0  |  v1.1.0  | Code used for the 2015 Run 5-6 electron-antineutrino appearance analysis.
             |          | This version was used for the first version of the tech note.
             |          | (2.5E20 POT, data not yet fit)
             |          | Details may be found in T2K-TN-252v1.
             |          |
2016_v1.0.0  |  v1.2.0  | Code used for the 2016 Run 1-6 / Run 1-7 electron-antineutrino appearance analysis.
             |          | Details may be found in T2K-TN-296.
             |          |
2016_v1.0.0  |  v1.2.0  | Code used for the 2016 Run 1-6 / Run 1-7 3-flavour muon-neutrino/muon-antineutrino disappearance analysis.
             |          | Details may be found in T2K-TN-302.
             |          |
2016_v1.0.0  |  v1.2.0  | Code used for the 2016 Run 1-6 / Run 1-7 3-flavour joint electron-neutrino/electron-antineutrino 
             |          | appearance and muon-neutrino/muon-antineutrino disappearance analysis. 
             |          | This analysis uses the standard 4 samples:
             |          | Single-ring e-like and single-ring mu-like samples, both in neutrino and antineutrino modes.
             |          | Details may be found in T2K-TN-266.
             |          |        
2016_v1.0.0  |  v1.2.0  | Code used for the 2016 Run 1-6 / Run 1-7 3-flavour joint electron-neutrino/electron-antineutrino 
             |          | appearance and muon-neutrino/muon-antineutrino disappearance analysis. 
             |          | This analysis incorporates a 5th SuperK (electron-neutrino CC1\pi^+ enhanced) sample 
             |          | into the VALOR T2K oscillation analysis.
             |          | Details may be found in T2K-TN-306.
             |          |        
2017_v1.0.0  |  v1.4.2  | Code used for the 2017 Run 1-8 3-flavour joint electron-neutrino/electron-antineutrino 
             |          | appearance and muon-neutrino/muon-antineutrino disappearance analysis (5 SuperK samples).
             |          | Details may be found in T2K-TN-327.
             |          |        
2018_v1.0.0  |  v1.5.0  | Code used for the 2018 Run 1-9c 3-flavour joint electron-neutrino/electron-antineutrino 
             |          | appearance and muon-neutrino/muon-antineutrino disappearance analysis (5 SuperK samples).
             |          | Also used for the 2018 Run 1-9c electron-antineutrino appearance analysis.
             |          | Details may be found in the 2018 statistical update to T2K-TN-327 and in T2K-TN-356.
             |          |        
2019_v1.0.0  |  v1.6.0  | Code used for the 2019 Run 1-9d 3-flavour joint electron-neutrino/electron-antineutrino 
             |          | appearance and muon-neutrino/muon-antineutrino disappearance analysis (5 SuperK samples).
             |          | Also used for the 2019 Run 1-9d electron-antineutrino appearance analysis.
             |          | Details may be found in T2K-TN-360 and T2K-TN-356.
             |          |       
?            |  ?       | Code used for the 2020 Run 1-10 3-flavour joint electron-neutrino/electron-antineutrino 
             |          | appearance and muon-neutrino/muon-antineutrino disappearance analysis (5 SuperK samples).
             |          | Details may be found in T2K-TN-394.
            
</pre>
