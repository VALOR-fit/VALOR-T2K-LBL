#!/bin/env python

import argparse
import ROOT
import sys

delim_list = lambda s: sorted(list(set(s.split(','))))
pair = lambda s: s.split(',') if len(s.split(',')) == 2 else sys.exit(s)
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--param0',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--param1',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--param2',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--param3',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--param4',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--param5',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--param6',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--param7',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--param8',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--param9',type=pair,default=None,help='Use comma separated list of postfit covariance and central value')
parser.add_argument('--increasenueflux',default=None,type=float,help='Increase nue flux by this amount. 0.02 corresponds to a 2%% increase')
parser.add_argument('--extrafluxfile',default=None,help='A file that has a flux matrix in (also with ND bins...) to be added to the input BANFF file values')
parser.add_argument('--infile', '-i',type=str,default='$T2KOSCANA/inputs/official/2018/banff/DataFit_Postfit_2018_final_v1_180907_sk_eb_valor_order_all_plus_scc.root',help='Input BANFF file')
parser.add_argument('--outfiletag','-o',type=str,default='$T2KOSCANA/inputs/hk/banff/banff.root',help='New location & tag')
args = parser.parse_args()

print parser
print args

#setup VALOR
#ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");

ROOT.gROOT.SetBatch(1)

valor = ROOT.valor
util  = valor.util
conv  = valor.conv

#get the original objects
banff_file = ROOT.TFile(args.infile,"read")
param_list = banff_file.Get("param_list")
postfit_cov = banff_file.Get("postfit_cov")
prefit_cov  = banff_file.Get("prefit_cov")
postfit_params = banff_file.Get("postfit_params")
prefit_params  = banff_file.Get("prefit_params")
banff_untuned_central_values = banff_file.Get("banff_untuned_central_values")
param_list_banff_xsec = banff_file.Get("param_list_banff_xsec")
banff_xsec_untuned_central_values = banff_file.Get("banff_xsec_untuned_central_values")

#get sizes
nrows = postfit_cov.GetNrows()
ncols = postfit_cov.GetNcols()
assert(nrows == ncols)
new_params = {}
inewparam = 0
if args.param0:
    new_params[inewparam] = {'param':0, 'error':float(args.param0[0]), 'central':float(args.param0[1])}
    inewparam += 1
if args.param1:
    new_params[inewparam] = {'param':1, 'error':float(args.param1[0]), 'central':float(args.param1[1])}
    inewparam += 1
if args.param2:
    new_params[inewparam] = {'param':2, 'error':float(args.param2[0]), 'central':float(args.param2[1])}
    inewparam += 1
if args.param3:
    new_params[inewparam] = {'param':3, 'error':float(args.param3[0]), 'central':float(args.param3[1])}
    inewparam += 1
if args.param4:
    new_params[inewparam] = {'param':4, 'error':float(args.param4[0]), 'central':float(args.param4[1])}
    inewparam += 1
if args.param5:
    new_params[inewparam] = {'param':5, 'error':float(args.param5[0]), 'central':float(args.param5[1])}
    inewparam += 1
if args.param6:
    new_params[inewparam] = {'param':6, 'error':float(args.param6[0]), 'central':float(args.param6[1])}
    inewparam += 1
if args.param7:
    new_params[inewparam] = {'param':7, 'error':float(args.param7[0]), 'central':float(args.param7[1])}
    inewparam += 1
if args.param8:
    new_params[inewparam] = {'param':8, 'error':float(args.param8[0]), 'central':float(args.param8[1])}
    inewparam += 1
if args.param9:
    new_params[inewparam] = {'param':9, 'error':float(args.param9[0]), 'central':float(args.param9[1])}
    inewparam += 1
num_new_params = nrows + len(new_params)
num_new_xsec_params = param_list_banff_xsec.GetEntriesFast() + len(new_params)

#get the additional flux matrix
if args.extrafluxfile:
    flux_file = ROOT.TFile(args.extrafluxfile, "read")
    flux_matrix = flux_file.Get("total_fdloc_cov")
    nrows_flux = flux_matrix.GetNrows()
    ncols_flux = flux_matrix.GetNcols()
    assert(nrows_flux == ncols_flux)
    assert(nrows_flux == 100)

#make new objects
output_param_list = ROOT.TObjArray()
output_postfit_cov = ROOT.TMatrixD(num_new_params, num_new_params)
output_prefit_cov = ROOT.TMatrixD(num_new_params, num_new_params)
output_prefit_params = ROOT.TVectorD(num_new_params)
output_postfit_params = ROOT.TVectorD(num_new_params)
output_param_list_banff_xsec = ROOT.TObjArray()
output_banff_untuned_central_values = ROOT.TVectorD(num_new_params)
output_banff_xsec_untuned_central_values = ROOT.TVectorD(num_new_xsec_params)

if args.extrafluxfile:
    f = open('difference_due_to_additional_flux_file.txt', 'w')

# Fill existing parameters
for inew in range(0,nrows):
    output_param_list.Add(ROOT.TObjString(param_list.At(inew).String().Data()))
    output_prefit_params[inew] = prefit_params[inew]
    output_postfit_params[inew] = postfit_params[inew]
    output_banff_untuned_central_values[inew] = banff_untuned_central_values[inew]
    for jnew in range(0,nrows):
        output_postfit_cov[inew][jnew] = postfit_cov(inew,jnew)
        output_prefit_cov [inew][jnew] = prefit_cov(inew,jnew)
        if args.extrafluxfile and inew < 50 and jnew < 50:
            if inew == jnew:
                x = ROOT.TMath.Sqrt(output_postfit_cov [inew][jnew])
            output_postfit_cov[inew][jnew] += flux_matrix(inew + 50, jnew + 50)
            output_prefit_cov [inew][jnew] += flux_matrix(inew + 50, jnew + 50)
            if inew == jnew:
                y = ROOT.TMath.Sqrt(output_postfit_cov [inew][jnew])
                f.write('Parameter ' + str(inew) + ' uncertainty has moved from ' + str(x) + ' to ' + str(y) + ' diff: ' + str(y-x) + '\n')
f.close()

#need to increase the nue flux in the numu beam (parameters 16-24)
#respect correlations too...
if args.increasenueflux:
    r = range(16,25)
    print r
    for inew in range(0,nrows):
        for jnew in range(0,nrows):
            if inew == jnew: #do diagonals at end
                continue
            if inew in r or jnew in r:
                oldcov = output_postfit_cov(inew,jnew)
                err_i  = ROOT.TMath.Sqrt(output_postfit_cov(inew,inew))
                err_j  = ROOT.TMath.Sqrt(output_postfit_cov(jnew,jnew))
                corr   = oldcov / (err_i * err_j)
                if inew in r:
                    err_i += args.increasenueflux
                if jnew in r:
                    err_j += args.increasenueflux
                newcov = corr * (err_i * err_j)
                output_postfit_cov[inew][jnew] = newcov
    for inew in r:
        err_i = ROOT.TMath.Sqrt(output_postfit_cov(inew,inew))
        print "Changing error in bin", inew, "from", err_i,
        err_i += args.increasenueflux
        print  "to", err_i
        newcov = err_i * err_i
        output_postfit_cov[inew][inew] = newcov

for inew in range(0,param_list_banff_xsec.GetEntriesFast()):        
    output_param_list_banff_xsec.Add(ROOT.TObjString(param_list_banff_xsec.At(inew).String().Data()))
    output_banff_xsec_untuned_central_values[inew] = banff_xsec_untuned_central_values[inew]

# Fill new parameters
for inew in range(0, len(new_params)):
    print util.AsStringCompact(conv.kNuisBANFFNew00 + inew)
    print ROOT.TObjString(util.AsStringCompact(conv.kNuisBANFFNew00 + inew))
    output_param_list.Add(ROOT.TObjString(util.AsStringCompact(conv.kNuisBANFFNew00 + inew)))
    print output_param_list[nrows + inew]
    output_prefit_params               [inew + nrows] = new_params[inew]['central']
    output_postfit_params              [inew + nrows] = new_params[inew]['central']
    output_banff_untuned_central_values[inew + nrows] = new_params[inew]['central']
    for jnew in range(0,len(new_params)):
        value = (ROOT.TMath.Power(new_params[inew]['error'], 2) if inew == jnew else 0)
        output_postfit_cov[inew + nrows][jnew + nrows] = value
        output_prefit_cov [inew + nrows][jnew + nrows] = value
    output_param_list_banff_xsec.Add(ROOT.TObjString(util.AsStringCompact(conv.kNuisBANFFNew00 + inew)))
    output_banff_xsec_untuned_central_values[inew + param_list_banff_xsec.GetEntriesFast()] = new_params[inew]['central']

#save
outfilename = args.outfiletag.split('.',1)[0]
for inew in range(0, len(new_params)):
    outfilename += '_%02d' % (new_params[inew]['param'])
if args.extrafluxfile:
    outfilename += '_extraflux'
if args.increasenueflux:
    outfilename += '_increasednueflux'
outfilename += '.root'
print outfilename
outfile = ROOT.TFile(outfilename,"recreate")
outfile.cd()
output_param_list.Write("param_list",ROOT.TObject.kSingleKey)
output_postfit_cov.Write("postfit_cov")
output_prefit_cov.Write("prefit_cov")
output_prefit_params.Write("prefit_params")
output_postfit_params.Write("postfit_params")
output_param_list_banff_xsec.Write("param_list_banff_xsec",ROOT.TObject.kSingleKey)
output_banff_untuned_central_values.Write("banff_untuned_central_values")
output_banff_xsec_untuned_central_values.Write("banff_xsec_untuned_central_values")
outfile.Close()
print 'Done'

#  LocalWords:  jnew
