#!/usr/bin/env python

'''

Script computing the maximum rate of fatal errors which could occur in the tracker due to
CIC faults, assuming that no such errors were observed during the total observation campaign

The estimated cross section limit is 1/(duration*fluence)

We estimate the limit for different LET threshold (E0=0.232*LET_thresh), and in order to determine
realistic cross section evolution we use a Weibull distribution for the integrated cross section, based on w and s computed on classic SEU data (see Compute_SEU_rates method)

Using a step function would be highly unrealistic and lead to artificially high values.

Values are computed for MPA and CBC mode, and are given in max number of upset per tracker per hour.

The particle rates used are identical to the Compute_SEU_rates method, and corresponds to a worst case scenario.

Options :
    1. -m/--mode       -> data mode (CBC or MPA)


'''

import sys, getopt
import os
import numpy as np
import scipy as sp
import scipy.stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
from array import array
from lmfit import Minimizer, Parameters, report_fit
from matplotlib.colors import LogNorm
from matplotlib import ticker
import matplotlib.colors as colors



# Here it's for the rate calculation
def func(x,a,b,c,d):
  #print ("In Function:", x,a,b,c,d )
  if c!=0:
    y = a*(1-np.exp(-((x - b)/c)**d))
  else:
    y = 0
  return y

# Calculation of the overall cross section @LHC once the device cross section is known
def SEU_xsec_lim(thresh,sigma,w,s):
  sigmas = [] # Device cross section
  probb = []  # LHC cross section

  for row in lineprob:
    stuff=row.split('\n')[0].split(',')
    if float(stuff[0])>thresh: # Xsec=0 in this case
        sigmas.append(func(float(stuff[0]), sigma,thresh,w,s))
    else:
        sigmas.append(0.)
    probb.append(float(stuff[1]))

  SEU_xs=0
  for i in range(len(sigmas)):
    if i <len(probb)-1:
      SEU_xs += (sigmas[i+1]-sigmas[i])*probb[i]*1e+08 # probb is the SEU prob in a sensitive volume of 1 cubic micrometer, the cross section is on 1 square micrometer, so 10-8 cm-2

  return SEU_xs

'''
The standard macro starts here
'''

try:
    opts, args = getopt.getopt(sys.argv[1:],"hm:",["mode="])
except getopt.GetoptError:
    outlog.write('UpsetRate.py -m <Mode>')
    sys.exit(2)

mode=''
verbose=0

for opt, arg in opts:
    if opt == '-h':
        outlog.write('UpsetRate.py -m <Mode>')
        sys.exit()
    elif opt in ("-m", "--mode"):
        mode=arg
        print(arg)

'''
Define the parameters here
'''

fluence=15000   # The average fluence of all the recorded SEU data (in part/s)
duration=13300  # The total SEU run duration (without error)
s       = 1.    # s Weibull parameter from data fit
w       = 10.    # w Weibull parameter from data fit


hrate=6.24*1e6
nmods=7608
if mode == 'MPA':
    hrate=3.85e7
    nmods=5592

'''
Open the input files

datf is the text file containing the SEU run summary
probf is the table giving the probability of an upset for a given Edep, knowing
that Edep=0.232*LET
'''

probf=open('SEUrate.txt','r')


lineprob=probf.readlines()

prob=[]

# Retrieve the data

for line in lineprob:
    prob.append(line.split(','))

t = np.arange(3000., 10000., 100.)
let = np.arange(0.1, 5., 0.1)

z=[]

for i in range(len(let)):
  #print(i)
  toto=np.zeros(len(t))
  for j in range(len(t)):
    #print(i,j,let[i],1./t[j])
    toto[j]=3600*hrate*2*nmods*SEU_xsec_lim(let[i]*0.232,1./(t[j]*fluence),w,s)
  z.append(toto)

fig,ax=plt.subplots(1,1)

cp = ax.contour(t, let, np.asarray(z),50, cmap='gist_rainbow')
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Max number of fatal upset per hour (full detector)')
plt.clabel(cp, inline=1, fontsize=8, colors='k')
ax.set_ylabel("LET (in MeV/$cm^2$/mg)")
ax.set_xlabel("Time without errors (in s)")
plt.grid(True, which="both", ls="--")
rate_at_5   = 3600*hrate*2*nmods*SEU_xsec_lim(5.*.232,1./(duration*fluence),w,s)
rate_at_0p1 = 3600*hrate*2*nmods*SEU_xsec_lim(0.1*.232,1./(duration*fluence),w,s)

#print(SEU_xsec_lim(5.*0.232,1./(duration*fluence),w,s))
print("")
print("Compute the maximum rate of fatal upset",mode,"case: ")
print("Considering that no fatal errors were observe during the total length of the tests",duration,"s at a fluence of",fluence,"particles/s")
print("")
print("-> If LET threshold is 0.1 MeV/$cm^2$/mg limit is",'{:.1f}'.format(rate_at_0p1),"error/h")
print("-> If LET threshold is 5.0 MeV/$cm^2$/mg limit is",'{:.1f}'.format(rate_at_5),"error/h")

#plt.show()
plt.savefig("Hard_error_limit"+mode+".png")
