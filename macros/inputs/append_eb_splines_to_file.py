#!/usr/bin/env python

import ROOT
import sys

try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")
  
ROOT.gROOT.SetBatch(1)

# File naming convention
standard_pattern = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/data/official/2020/xsec/splines2020v5/spline2019splinesfitqunPre[NDIM]_t2ksk19b.fqv4r0.[CURRENT].[TYPE]_[SELECTION]selec[HASEB].root"
eb_pattern       = "/hepstore/fbench/software/XSecResponse/eb_splines[NDIM]_[TYPE]_[CURRENT]_[SELECTION]selec_2020-04-01.root"

def get_file_suffix(sample):
    if "1Rmu" in sample:
        return ""
    else:
        return "erectheta"


def main():
    beams           = ["fhc", "rhc"]
    osc_channels    = ["nuebar_x_nuebar", "nue_x_nue", "numubar_x_nuebar", "numubar_x_numubar", "numu_x_nue", "numu_x_numu"]
    selections      = ["numu", "nue", "cc1pi"]
    is_selection_2d = {"numu" : False, "nue" : True, "cc1pi" : True}

    for beam in beams:
        for channel in osc_channels:
            for selection in selections:
                fname_spline = standard_pattern
                fname_eb     = eb_pattern

                fname_spline_2d = ""
                fname_eb_2d     = ""
                if is_selection_2d[selection]:
                    fname_spline_2d = "erectheta"
                    fname_eb_2d     = "_2d"

                fname_spline = fname_spline.replace("[CURRENT]"  , beam)
                fname_spline = fname_spline.replace("[TYPE]"     , channel)
                fname_spline = fname_spline.replace("[SELECTION]", selection)
                fname_spline = fname_spline.replace("[NDIM]"     , fname_spline_2d)

                fname_out    = fname_spline.replace("[HASEB]", "_plus_Eb")
                fname_spline = fname_spline.replace("[HASEB]", "")

                fname_eb = fname_eb.replace("[CURRENT]"  , beam)
                fname_eb = fname_eb.replace("[TYPE]"     , channel)
                fname_eb = fname_eb.replace("[SELECTION]", selection)
                fname_eb = fname_eb.replace("[NDIM]"     , fname_eb_2d)

                fspline = ROOT.TFile(fname_spline, "READ")
                feb     = ROOT.TFile(fname_eb    , "READ")
                fout    = ROOT.TFile(fname_out   , "RECREATE")

                if not fspline.IsOpen():
                    print "ERROR: Could not find the file [{}]".format(fname_spline)
                    sys.exit()
                if not feb.IsOpen():
                    print "ERROR: Could not find the file [{}]".format(fname_eb)
                    sys.exit()
                if not fout.IsOpen():
                    print "ERROR: Could not create the file [{}]".format(fname_out)
                    sys.exit()

                for file in [fspline, feb]:
                    #print "\n\n\n\n"
                    if file.GetListOfKeys().GetEntries() == 0:
                        print "ERROR: The input file [{}] is empty!".format(file.GetName())
                        sys.exit()
                    for key in file.GetListOfKeys():
                        #print "Copying", key.GetName(), "from", file.GetName(), "to", fout.GetName()
                        fout.WriteObject(key.ReadObj(), key.GetName())

                print "\n\nCopied", fspline.GetListOfKeys().GetEntries(), "splines from", fspline.GetName(), "to", fout.GetName()
                print "Copied", feb    .GetListOfKeys().GetEntries(), "splines from", feb    .GetName(), "to", fout.GetName()

                fout.Close()
                feb.Close()
                fspline.Close()

    print "Done!"

  
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="%prog [options]..",description="Append Eb splines to the standard spline files")
    parser.add_argument('-t', '--tag', default='spline2019v1', help="Spline filename prefix.")
    args = parser.parse_args()
    main()
 
