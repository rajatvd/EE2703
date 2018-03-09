
from pylab import *
import sys
import argparse

ARGS = {"n":100,
       "M":5,
       "nk":500,
       "u0":5,
       "p":0.25}

def addargs(args,parser):
    for arg in args:
        parser.add_argument('-{}'.format(arg),default=args[arg],type=type(args[arg]))
        
ap = argparse.ArgumentParser()
addargs(ARGS,ap)
args = ap.parse_args()

for arg in ARGS:
    exec("{} = args.{}".format(arg,arg))
    
print(n)