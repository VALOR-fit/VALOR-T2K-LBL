#!/usr/bin/env python

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)

samples=["1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe"]

def convert_hybrid(outfile, datafile, fakefile, ntoys, data_samples=[]):
    data = ROOT.TFile(datafile, "READ")
    hybrid = ROOT.TFile(outfile, "RECREATE")
    fake = ROOT.TFile(fakefile, "READ")
    print "Reading fake data toys from", fakefile, ", output hybrid to", outfile

    # Read the data spectra that we're interested in
    data_hists = {}
    data.cd()
    for sample in data_samples:
        data_hists[sample] = data.Get("mock_data_" + sample + "_0")

    for itoy in range (0, ntoys):
        if itoy % 1000 == 0:
            print "Processing toy", str(itoy)
        # Copy over the data and fake data spectra to the hybrid file
        for sample in samples:
            hname_suffix = "_data_" + sample + "_" + str(itoy)
            if sample in data_samples:
                # Replace the fake data spectrum with the data spectrum
                hybrid.cd()
                data_hists[sample].Write("mock" + hname_suffix)
                data_hists[sample].Write("exp" + hname_suffix)
            else:
                # Copy the fake data spectrum to the output file
                fake.cd()
                hfake_mock = fake.Get("mock" + hname_suffix)
                hfake_exp = fake.Get("exp" + hname_suffix)
                hybrid.cd()
                hfake_mock.Write("mock" + hname_suffix)
                hfake_exp.Write("exp" + hname_suffix)
                del hfake_mock
                del hfake_exp

    # Copy over all other objects from the fake data file to the hybrid file
    for path, object in getall(fake):
        if not("mock_data_" in path or "exp_data_" in path):
            hybrid.cd()
            object.Write(path.replace("/","",1))

    data.Close()
    fake.Close()
    hybrid.Close()


def getall(d, basepath="/"):
    "Generator function to recurse into a ROOT file/dir and yield (path, obj) pairs"
    d.cd()
    for key in d.GetListOfKeys():
        kname = key.GetName()
        if key.IsFolder():
            for i in getall(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield basepath+kname, d.Get(kname)


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Convert and already existing set of fake data toys into hybrid toys")
    parser.add_option("-o", "--output-file", dest="outfile", help='Output file name for hybrid toys')
    parser.add_option("-d", "--data-file", dest="datafile", help='Input file containing data spectra')
    parser.add_option("-f", "--fake-data", dest="fakefile", help='Full path of fake data toys')
    parser.add_option('-n', "--ntoys", dest='ntoys', default=1, type=int, help='Number of fake data toys in the input file')
    
    v = optparse.Values()
    parser.parse_args(values = v)
    o, other = parser.parse_args()

    convert_hybrid(o.outfile, o.datafile, o.fakefile, o.ntoys, ["1Rmu","1Re","RHC1Rmu","MultiRe"])

    if other:
        print "Warning - ignored arguments: ", other
