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

def main(infilename, infilename_data, ntoys, nthrows, tag, itoy_start, analytic, reweighting, chisq_min):
    # Note: in this macro we deal entirely with calculating an unweighted p-value for each toy. The weights should be applied at the next step, i.e. when calculating the (weighted) mean p-value
    # Set up outputs
    Nevents_max = 100
    outname = tag + ".root"
    fout = ROOT.TFile(outname, "recreate")
    tree_out = ROOT.TTree("output_tree","output_tree")
    pval_left         = array('f', [-1.])   # For why an array is used see https://root.cern.ch/how/how-write-ttree-python
    pval_right        = array('f', [-1.])
    Nev               = array('f', [-1.])
    hist_nevents_unweighted_itoy = ROOT.TH1D("Nev_fluct_itoy", "Nev_fluct_itoy", Nevents_max, 0, Nevents_max)
    Nev_fluct_hist    = ROOT.TH1D("Nev_fluct", "Nev_fluct", Nevents_max, 0, Nevents_max)     # store all the throws in this hist rather than a vector, due to storage space

    tree_out.Branch("pvalue_left"       , pval_left        , "pvalue_left/F" )
    tree_out.Branch("pvalue_right"      , pval_right       , "pvalue_right/F")
    tree_out.Branch("Nevents"           , Nev              , "Nevents/F"     )
    tree_out.Branch("hist_nevents_fluct", hist_nevents_unweighted_itoy       )

    # Read inputs
    fin = ROOT.TFile(infilename, "read")
    if not fin:
        print "Cannot open fake-data file", infilename, ". Exiting."
        util.Die(-1)
    fin_data = ROOT.TFile(infilename_data, "read")
    if not fin_data:
        print "Cannot open data file", infilename_data, ". Exiting."
        util.Die(-1)
    hdata_RHC1Re  = fin_data.Get("mock_data_RHC1Re_0")
    Nevents_data  = hdata_RHC1Re.Integral()
    print "start toy =", itoy_start
    print "ntoys =", ntoys
    for itoy in range(itoy_start, itoy_start + ntoys):
        if itoy % 1000 == 0:
            print "itoy =", itoy, "of", ntoys + itoy_start

        list_nevents_unweighted_itoy = []

        hname = "exp_data_RHC1Re_"  + str(itoy)
        hErecoToy_RHC1Re = fin.Get(hname)
        if not hErecoToy_RHC1Re:
            print "WARNING: Histogram", hname, "not found in input file. Skipping this toy."
            continue
        Nevents_nuebar = hErecoToy_RHC1Re.Integral()
        Nev[0] = Nevents_nuebar

        if analytic == 0:
            # Fluctuate the nuebar event rate
            for ithrow in range(0, nthrows):
                Nevents_fluct = rnd.RndToyMCstat().Poisson(Nevents_nuebar)
                list_nevents_unweighted_itoy.append(Nevents_fluct)
                hist_nevents_unweighted_itoy.Fill(Nevents_fluct)
                Nev_fluct_hist.Fill(Nevents_fluct)

            # Calculate p-values for this fake-data toy
            # Note: the weight is not applied when calculating rate-only p-values, as these p-vals are claculated from
            # stat throws of the RHC1Re sample, but the weight is only dependant on the control samples, so would cancel out
            pval_left [0] = p_value_from_list(list_nevents_unweighted_itoy, Nevents_data, "left")
            pval_right[0] = p_value_from_list(list_nevents_unweighted_itoy, Nevents_data, "right")
        else:
            pval_left [0] = p_value_analytic(Nevents_nuebar, math.floor(Nevents_data), "left")
            pval_right[0] = p_value_analytic(Nevents_nuebar, math.ceil(Nevents_data), "right")

        tree_out.Fill()
        del list_nevents_unweighted_itoy
        hist_nevents_unweighted_itoy.Reset()

    # Write out
    fout.cd()
    Nev_fluct_hist.Write("Nevents_nuebar_fluct")
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

def p_value_from_list(list, data, tail):
    pval = 0
    for Nev in list:
        if tail == "right":
            if Nev >= data:
                pval += 1
        elif tail == "left":
            if Nev <= data:
                pval += 1
        else:
            print "p_value: Unrecognised argument", tail, "exiting."
            util.Die(-1)
    return float(pval) / float (len(list))

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
    parser.add_option('--reweighting', dest='reweighting', action='store_true', default="False", help='Include this flag to explicitly use control sample weights in p-value calculations. Only use if fake data toys were not generated by MC sampling')
    parser.add_option('--chisq-min', dest='chisq_min', default=0, type=float, help='Minimum control sample chisq')


    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.infilename, o.infilename_data, o.ntoys, o.nthrows, o.tag, o.itoy_start, o.analytic, o.reweighting, o.chisq_min)

    if other:
        print "Warning - ignored arguments: ",other
                                                   
