#!/bin/env python2

from argparse import ArgumentParser
from shutil import copyfile
from os.path import expandvars

parser = ArgumentParser()
parser.add_argument('-c','--config-file',required=True,type=str)
parser.add_argument('-d','--date',required=True,type=str)
args = parser.parse_args()

with open(expandvars('$T2KOSCANA/config/' + args.config_file)) as f:
    foundbanff = False
    for l in f.readlines():
        l = l.strip()
        if not l.startswith('<filename>'):
            continue
        l = l.split('filename>')[1].strip()
        if not l:
            continue
        l = l.split('<')[0].strip()
        if '/skdetfsi/' in l or '/banff/' in l:
            if '/banff/' in l:
                if foundbanff:
                    continue
                foundbanff = True
            orig = expandvars('$T2KOSCANA/inputs/' + l)
            new  = expandvars('$HOME/public_html/oa/' + args.date + '/$systtag/' + l.split('/')[-1])
            print expandvars('$T2KOSCANA/inputs/' + l), new
            copyfile(orig, new)
