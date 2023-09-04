#!/usr/bin/env python

import ROOT
import sys

try:
    ROOT.valor

except:
    ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

  
ROOT.gROOT.SetBatch(1)

# File naming convention
output_pattern   = "$RUN/appended/spline2021splinesfitqunPre[NDIM]_t2ksk19b.fqv4r0.21bv2.[CURRENT].[TYPE]_[SELECTION]selec_coulomb_plus_Eb.root"
standard_pattern = "$RUN/raw_splines/xsec_splines[NDIM]_coulomb_[TYPE]_[CURRENT]_[SELECTION]selec_2023_05_08.root"
eb_pattern       = "$RUN/raw_splines/eb_splines[NDIM]_coulomb_[TYPE]_[CURRENT]_[SELECTION]selec_2023_05_08.root"

def main():
   # beams           = ["fhc", "rhc"]
    beams           = ["fhc"]
   # beams           = ["rhc"]
    osc_channels    = ["nuebar_x_nuebar", "nue_x_nue", "numubar_x_nuebar", "numubar_x_numubar", "numu_x_nue", "numu_x_numu"]
   # selections      = ["numu", "nue"]
    selections      = ["numu", "nue", "cc1pi", "numucc1pi"]
    is_selection_2d = {"numu" : False, "nue" : True, "cc1pi" : True, "numucc1pi": False}

    for beam in beams:
        for channel in osc_channels:
            for selection in selections:
                if beam =="rhc" and (selection =="numucc1pi"): #skipping rch for multirmu as this is FHC sample 
                   print "skipping ", selection, " in channel ", channel
                   continue
                fname_spline = standard_pattern
                fname_eb     = eb_pattern
                fname_out    = output_pattern

                fname_valor_2d        = ""
                fname_xsecresponse_2d = ""
                if is_selection_2d[selection]:
                    fname_valor_2d        = "erectheta"
                    fname_xsecresponse_2d = "_2d"

                fname_spline = fname_spline.replace("[CURRENT]"  , beam)
                fname_spline = fname_spline.replace("[TYPE]"     , channel)
                fname_spline = fname_spline.replace("[SELECTION]", selection)
                fname_spline = fname_spline.replace("[NDIM]"     , fname_xsecresponse_2d)

                fname_eb = fname_eb.replace("[CURRENT]"  , beam)
                fname_eb = fname_eb.replace("[TYPE]"     , channel)
                fname_eb = fname_eb.replace("[SELECTION]", selection)
                fname_eb = fname_eb.replace("[NDIM]"     , fname_xsecresponse_2d)

                fname_out = fname_out.replace("[CURRENT]"  , beam)
                fname_out = fname_out.replace("[TYPE]"     , channel)
                fname_out = fname_out.replace("[SELECTION]", selection)
                fname_out = fname_out.replace("[NDIM]"     , fname_valor_2d)

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
 
