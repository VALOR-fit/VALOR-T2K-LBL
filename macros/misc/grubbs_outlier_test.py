import ROOT
import numpy

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

def main(str_list, alpha):
    # See https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h1.htm for explanation of this test

    # Convert input list of strings into list of floats and sort into ascending order
    values = [float(item) for item in str_list]
    values.sort()
    print "input list =", values

    # Test if min value is an outlier
    test_val = min(values)
    grubbs_test(test_val, values, alpha)

    # Test if max value is an outlier
    test_val = max(values)
    grubbs_test(test_val, values, alpha)

    print "Done!"


def grubbs_test(test_value, values, alpha):
    # construct a sample that does not include the tested value
    sample = values
    sample.remove(test_value)
    print "For tested value", test_value, "and sample", sample, ":"

    # Calculate Grubb's test statistic
    statistic = abs(test_value - numpy.mean(sample)) / numpy.std(sample, ddof=1)

    critval = get_2side_t_critval(sample, alpha)

    # Calculate rejection threshold and compare to test statistic
    threshold = ( float(len(sample) - 1.) / ROOT.TMath.Sqrt(float(len(sample))) ) * ROOT.TMath.Sqrt(critval * critval / ( float(len(sample)) - 2. + critval * critval ))
    print "H0:  There are no outliers in the data"
    print "Ha:  The value", test_value, "is an outlier"
    print "Grubb's test statistic, G =", statistic
    print "Significance level =", alpha
    print "Critical value for a two-tailed test =", critval
    print "Critical region:  Reject H0 if G >", critval
    if statistic > critval:
        print "Outcome: H0 rejected in favour of Ha - this value is an outlier at this significance"
    else:
        print "Outcome: H0 not rejected - this value is not outlier at this significance"
    print "\n\n"


def get_2side_t_critval(sample, alpha):
    prob = 1. - ( alpha / float(2. * len(sample)) )
    dof = len(sample) - 2
    return ROOT.TMath.StudentQuantile(prob, dof)


def get_comma_separated_args(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Perform Grubb's test for outliers on an input list of values")
    parser.add_option('--values', dest='values', type='string', action='callback', callback=get_comma_separated_args, help='List of values to test for outliers')
    parser.add_option('--significance', dest='alpha', default=0.05, type=float, help='Significance level used to reject hypothesis')

    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.values, o.alpha)

    if other:
        print "Warning - ignored arguments: ",other
                                                   
