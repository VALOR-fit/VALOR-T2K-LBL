#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging
import glob # file matching
import tarfile

#loop through directories
rootdir = "/Users/thomasdealtry/T2K/ralrepo/t2k/t2kana/nuosc/trunk/2012results/"
print rootdir
for root, subFolders, files in os.walk(rootdir):
    print "root"
    print root
    print "subFolders"
    print subFolders
    #print files

    #find all pbs, pbs.err, pbs.out files
    allpbsfiles = glob.glob(os.path.join(root,"*.pbs*"))

    #loop over file tags
    while allpbsfiles != []:
        print "allpbsfiles"
        print allpbsfiles
    
        #find the tag
        tarname = ""
        firstfilename = allpbsfiles[0].split('/')        
        firstname = firstfilename[-1].split('.pbs')
        tarname = firstname[0].split('-')
        print "tarname"
        print tarname
    
        #make a list with just these files
        thesepbsfiles = glob.glob(os.path.join(root,tarname[0]+"*.pbs*"))
        print "thesepbsfiles"
        print thesepbsfiles

        #change directory
        thisdir = thesepbsfiles[0].rsplit('/',1)[0]
        os.chdir(thisdir)
        
        #tar 
        print "tar.gz filename"
        print os.path.join(root,tarname[0]+".tar.gz")
        tar = tarfile.open(os.path.join(root,tarname[0]+".tar.gz"), "w:gz")
        for thisname in thesepbsfiles:
            thisfilename = thisname.split('/')[-1]
            tar.add(thisfilename)
        tar.close()

        #remove the files we've just tarred
        tar = tarfile.open(os.path.join(root,tarname[0]+".tar.gz"), "r:gz")
        print "names in tar file"
        print tar.getnames()
        for thisname in tar.getnames():
            os.remove(thisname)
        tar.close()

        #find all pbs, pbs.err, pbs.out files again
        allpbsfiles = glob.glob(os.path.join(root,"*.pbs*"))

