import ROOT
import math
from array import array

ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");

ROOT.gROOT.SetBatch(1)

valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
t2kopt  = valor.T2KOpt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

def main(infilename, infilename_data, ntoys, nthrows, tag, itoy_start, analytic):
    # Note: in this macro we deal entirely with calculating an unweighted p-value for each toy. The weights should be applied at the next step, i.e. when calculating the (weighted) mean p-value
    # Set up outputs
    Nevents_max = 100
    outname = tag + ".root"
    fout = ROOT.TFile(outname, "recreate")
    hist_nevents = ROOT.TH1D("Nevents", "Nevents", Nevents_max, 0, Nevents_max)
    tree_out = ROOT.TTree("output_tree","output_tree")
    pval_left         = array('f', [-1.])   # For why an array is used see https://root.cern.ch/how/how-write-ttree-python
    pval_right        = array('f', [-1.])

    tree_out.Branch("pvalue_left"       , pval_left        , "pvalue_left/F"         )
    tree_out.Branch("pvalue_right"      , pval_right       , "pvalue_right/F"        )

    # Read inputs
    fin = ROOT.TFile(infilename, "read")
    if not fin:
        print "Cannot open fake-data file", infilename, ". Exiting."
        util.Die(-1)
    fin_data = ROOT.TFile(infilename_data, "read")
    if not fin_data:
        print "Cannot open data file", infilename_data, ". Exiting."
        util.Die(-1)
    Nevents_data_hist = fin_data.Get("mock_data_RHC1Re_0")
    Nevents_data = Nevents_data_hist.Integral()
    for itoy in range(itoy_start, itoy_start + ntoys):
        if itoy % 1000 == 0:
            print "itoy =", itoy - itoy_start, "of", ntoys

        hist_nevents_unweighted_itoy = ROOT.TH1D("Nevents_unweighted_" + str(itoy), "Nevents_unweighted_" + str(itoy), Nevents_max, 0, Nevents_max)

        hErecoToy_RHC1Re = fin.Get("mock_data_RHC1Re_"  + str(itoy))
        Nevents_nuebar = hErecoToy_RHC1Re.Integral()

        if analytic == 0:
            # Fluctuate the nuebar event rate
            for ithrow in range(0, nthrows):
                Nevents_fluct = rnd.RndToyMCstat().Poisson(Nevents_nuebar)
                hist_nevents_unweighted_itoy.Fill(Nevents_fluct)

            # Calculate p-values for this fake-data toy
            pval_left [0] = p_value(hist_nevents_unweighted_itoy, math.floor(Nevents_data), "left" )
            pval_right[0] = p_value(hist_nevents_unweighted_itoy, math.ceil(Nevents_data) , "right")
        else:
            pval_left [0] = p_value_analytic(Nevents_nuebar, math.floor(Nevents_data), "left")
            pval_right[0] = p_value_analytic(Nevents_nuebar, math.ceil(Nevents_data), "right")

        tree_out.Fill()
        del hist_nevents_unweighted_itoy

    # Write out
    fout.cd()
    hist_nevents.Write("Nevents_nuebar_unweighted")
    tree_out.Write("fake-data_tree")
    fout.Close()
    fin_data.Close()
    print "Done!"


def p_value(distribution, data, tail):
    ntot = distribution.Integral()
    ibin_data = distribution.FindBin(data)
    pval = -999
    if tail == "right":
        pval = distribution.Integral(ibin_data, distribution.GetNbinsX())
    elif tail == "left":
        pval = distribution.Integral(1, ibin_data)
    else:
        print "p_value: Unrecognised argument", tail, "exiting."
        util.Die(-1)
    return pval/ntot

def p_value_analytic(mean, data, tail):
    if tail == "right":
        return 1 - ROOT.Math.inc_gamma_c(data, mean)
    elif tail == "left":
        return ROOT.Math.inc_gamma_c(data + 1, mean)
    else:
        print "p_value: Unrecognised argument", tail, "exiting."
        util.Die(-1)


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Perform Poisson throws (default 10k) on the number of nuebar events for each of the fake-data toys in input file and store the weighted and unweighted flucts")
    parser.add_option('--infilename', dest='infilename', default="", help='Root file containing the fake-data toys to be fluctuated')
    parser.add_option('--infilename-data', dest='infilename_data', default="", help='Root file containing the data or the Asimov, in order to compute p-values')
    parser.add_option('--ntoys', dest='ntoys', default=100000, type=int, help='Number of fake-data toys in the input file')
    parser.add_option('--nthrows', dest='nthrows', default=10000, type=int, help='Number of Poisson throws')
    parser.add_option('--itoy-start', dest='itoy_start', default=0, type=int, help='Start reading the input file at this toy number')
    parser.add_option('--outfilename', dest='tag', default="nuebarapp_deltachisq", help='Output file name')
    parser.add_option('--analytic', dest='analytic', default=0, type=int, help='Analytically calculate the p-values instead of throwing many Poisson fluctuations')


    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.infilename, o.infilename_data, o.ntoys, o.nthrows, o.tag, o.itoy_start, o.analytic)

    if other:
        print "Warning - ignored arguments: ",other
                                                   
