#! /usr/bin/env python

"""
Used to get a good look at the MINUIT fit results

Takes a list of .out files from a production run that contain the standard MINUIT outputs.
Performs a grep, making a temporary file "tempfile.txt" in the working directory.
Runs over tempfile.txt to extract the important quantities:
  chisq, #calls, status, edm, strategy, err, fit parameter: final value, err, step size, first derivative
and put them into a TTree in file "minuit.root"
Run with:
 ./minuit_to_tree "/path/to/batch/outputs/*.out"
 (need quotes otherwise will only grep the first file the * returns)
 
Currently only works for 3 flavour osc. (i.e. fits with 6 parameters). Let me know if you want it more flexible & I can put something together relatively easily
The python is a bit of a hack too. Still learning...

T.Dealtry
"""


import shutil
import sys
import optparse
import os
import glob
import re
import math # cos, log, pi etc
import logging # print with logging.debug/info/warning/error/critical('message')
from array import array
import ROOT

#remove the temporary file, then run the grep on the input variables
#makes a file with lines that contain either:
#FCN (chisq), STATUS, CALLS
#EDM, STRATEGY, ERR MATRIX STATUS
#PARAM: value, err, step size, first derivative
#refit
#com = 'grep \"FCN=.*\..* FROM MIGRAD\|EDM\|1  sinsq2theta23\|2  sinsq2theta13\|3  sinsq2theta12\|4  delta_cp\|5  |dm23sq|\|6  dm12sq\|refit\" '
com = 'grep \"FCN=.*\..* FROM MIGRAD\|EDM\|^   [0-9]  \|^  [0-9][0-9]  \|refit\" --with-filename '
com = com + sys.argv[1]
com = com +  ' >> tempfile.txt'
print com
os.system('rm tempfile.txt')
os.system(com)

    
#open the output root file & set up the tree
tf = ROOT.TFile("minuit.root","recreate")
t = ROOT.TTree('fitResults','minuit fit results')

fcn = array('f',[-1.])
calls = array('i',[0])
status = array('i',[-1])
edm = array('f',[-1.])
strategy = array('i',[-1])
err = array('f',[-5.])
arraysize = 4
sinsq2theta23 = array('f',arraysize*[0.])
sinsq2theta13 = array('f',arraysize*[0.])
sinsq2theta12 = array('f',arraysize*[0.])
delta_cp = array('f',arraysize*[0.])
dm23sq = array('f',arraysize*[0.])
dm12sq = array('f',arraysize*[0.])

t.Branch('chisq',fcn,'chisq/F')
t.Branch('calls',calls,'calls/I')
t.Branch('status',status,'status/I')
t.Branch('edm',edm,'edm/F')
t.Branch('strategy',strategy,'strategy/I')
t.Branch('err',err,'err/F')
#arrays of size 4 to hold parameter:
# value, error, step size, first derivative
t.Branch('sinsq2theta23',sinsq2theta23,'sinsq2theta23[4]/F')
t.Branch('sinsq2theta13',sinsq2theta13,'sinsq2theta13[4]/F')
t.Branch('sinsq2theta12',sinsq2theta12,'sinsq2theta12[4]/F')
t.Branch('delta_cp',delta_cp,'delta_cp[4]/F')
t.Branch('dm23sq',dm23sq,'dm23sq[4]/F')
t.Branch('dm12sq',dm12sq,'dm12sq[4]/F')


def GetParamValues(line):
    """
    This function decodes a "parameter" line, returns an array of size 4
    """
    param = array('f',4*[0.])
    temp = line.split(None,6)
    param[0] = float(temp[3])
    if 'fixed' in temp[4]:
        param[1] = 0.
        param[2] = 0.
        param[3] = 0.
    else:
        param[1] = float(temp[4])
        param[2] = float(temp[5].strip('*'))
        if 'at limit **' in temp[6]:
            param[3] = -99.
        else:
            param[3] = float(temp[6])
    return param


#setup a few counters
nconverged = 0
nfailed = 0
ncalllimit = 0
nunknown = 0
#and a couple of switches used to navigate the file
skipline = True # look for first 'fcn' line
filltree = False # don't fill until we've got the values first.

#open the file & iterate over lines
f = open('tempfile.txt', 'r')
for line in f:
    #this is always the starting line of a block we want to save.
    if 'FCN' in line:
        #fill the tree, unless 1st line of file or refit line found above.
        if filltree:
            t.Fill()
        #reset some variables
        fcn[0] = -1.
        calls[0] = 0
        status[0] = -1
        #split up the line
        temp = line.split()
        #get the fcn
        fcntemp = temp[1].split("=")
        fcn[0] = float(fcntemp[1])
        #get the status
        if 'CONVERGED' in temp[4]:
            nconverged = nconverged + 1
            status[0] = 0
        elif 'FAILED' in temp[4]:
            nfailed = nfailed + 1
            status[0] = 1
            #print "failed", line
        elif 'CALL' in temp[4]:
            ncalllimit = ncalllimit + 1
            status[0] = 2
            print "call limit reached", line
        else:
            nunknown = nunknown + 1
            status[0] = 3
            print "unknown", line
        #get the number of calls
        calls[0] = int(temp[len(temp)-4])
        skipline = False
        continue
    #if we've got a refit, we don't want to save the previous block into the tree
    elif 'refit' in line:
        print line
        filltree = False
    #we'll keep skipping until we find another 'fcn' starting line
    if skipline:
        continue
    if 'EDM' in line:
        #reset some variables
        edm[0] = -1.
        strategy[0] = -1
        err[0] = -5.
        #split up the line
        temp = line.split(None, 6)
        #get the edm
        edmtemp = temp[1].split('=')
        edm[0] = float(edmtemp[1])
        #ge thte strategy
        strategy[0] = int(temp[3])
        #get the err matrix uncertainty
        if 'ACCURATE' in temp[6]:
            err[0] = 0.
        elif 'NOT POS-DEF' in temp[6]:
            err[0] = -1.
        elif 'APPROXIMATE' in temp [6]:
            err[0] = -2.
        else:
            errtemp = temp[6].split()
            err[0] = float(errtemp[1])
    elif 'sinsq2theta23' in line:
        holder = GetParamValues(line)
        #sinsq2theta23 = GetParamValues(line) doesn't fill the tree correctly.
        # no idea why... so need to do the following
        for i in xrange(4):
            sinsq2theta23[i] = holder[i]
    elif 'sinsq2theta13' in line:
        holder = GetParamValues(line)
        for i in xrange(4):
            sinsq2theta13[i] = holder[i]
    elif 'sinsq2theta12' in line:
        holder = GetParamValues(line)
        for i in xrange(4):
            sinsq2theta12[i] = holder[i]
    elif 'delta_cp' in line:
        holder = GetParamValues(line)
        for i in xrange(4):
            delta_cp[i] = holder[i]
    elif 'dm23sq' in line:
        holder = GetParamValues(line)
        for i in xrange(4):
            dm23sq[i] = holder[i]
    elif 'dm12sq' in line:
        holder = GetParamValues(line)
        for i in xrange(4):
            dm12sq[i] = holder[i]
        #dm12sq is the last in the block, so we want to skip forward to the next 'fcn' to fill the tree
        filltree = True
        skipline = True
        
#remember to fill the tree when we've reached the last line (there is no next 'fcn' starting line)
t.Fill()

#finish off by printing the counters & closing all the files.
print "converged:failed:callLimitReached:unknown", nconverged, nfailed, ncalllimit, nunknown
f.close()
#os.system('rm tempfile.txt')

tf.Write()
tf.Close()
