#!/bin/env python

import argparse
import os

parser = argparse.ArgumentParser(description='Run many WCSim jobs with different options. Use , to delimit multiple options', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#options about how to run this script
parser.add_argument('--onlycreatefiles', action='store_true', help="Do a test run where you create all the files, but don't submit?")
parser.add_argument('--notifyuseremail', type=str, default=None, help='Specify this to get email notifications about jobs')
#job setup options
parser.add_argument('--executable','-e', required=True, type=str, help='Executable')
parser.add_argument('--arguments', '-a', required=True, type=str, help='Arguments')
parser.add_argument('--logfile', '-l', required=True, type=str, help='Logfile')
parser.add_argument('--escape-quotes', default=False, action='store_true', help='Escape the " in the argument')
#parse the arguments
args = parser.parse_args()

fcondor = open(args.logfile + '.jdl', 'w')

condor = ''
if args.notifyuseremail:
    condor += '' \
        'notify_user    = ' + args.notifyuseremail + '\n' \
        'notification   = Always \n'
arguments = args.arguments
if args.escape_quotes:
    arguments = arguments.replace('"', '\\"')
condor += '' \
    'executable     = ' + args.executable + ' \n' \
    'universe       = vanilla \n' \
    'arguments      = ' + arguments + ' \n' \
    'output         = ' + args.logfile + '.out \n' \
    'error          = ' + args.logfile + '.err \n' \
    'log            = ' + args.logfile + '.log \n' \
    'request_memory = 1000 \n' \
    'getenv         = True \n' \
    'queue 1 \n'
fcondor.write(condor)
fcondor.close()
command = 'condor_submit ' + args.logfile + '.jdl'

#run the job
print command
if not args.onlycreatefiles:
    os.system(command)
else:
    print 'test run; not actually submitting'

#if __name__ == "__main__":
#    main()
        
