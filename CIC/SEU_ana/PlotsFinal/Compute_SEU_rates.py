#!/usr/bin/env python

'''

 Main script to compute the rates, according to the method
 described in the CIC SEU note:

 https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared%20Documents/Simulation%20studies/CIC2_SEU.pdf

The script retrieves the info, compute the cross sections for different LETs points and do the fit to the Weibull curve and finally computes the overall cross section per charged hadrons, and corresponding upset rate.

Errors are also computed, 95% Poisson CL for the number of errors and 5% for the ion fluxes

Options :
    1. -f/--file       -> Text file containing the test results for the different runs
                         (ordered in ascending LETs)
    2. -m/--mode       -> data mode (CBC or MPA)
    3. -t/--type       -> data path studied (L1 or Trigger)

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


# Function used for the fit (The Weibull curve)
def fcn2min(params, x, data):
    
    amp = params['sigma0']
    cen = params['E0']
    wid = params['W']
    d   = params['s']
    model = amp * (1-np.exp(-((x-cen)/wid)**d))
    return model - data
    
# Here it's for the rate calculation
def func(x,a,b,c,d):
  #print ("In Function:", x,a,b,c,d )
  if c!=0:
    y = a*(1-np.exp(-((x - b)/c)**d))
  else:
    y = 0
  #print (y)
  return y

# Calculation of the overall cross section @ LHC once the device cross section is known
def SEU_xsec(result):
  sigmas = [] # Device cross section
  probb = []  # LHC cross section

  for row in lineprob:
    stuff=row.split('\n')[0].split(',')
    if float(stuff[0])>result.params['E0'].value: # Xsec=0 in this case
        sigmas.append(func(float(stuff[0]), result.params['sigma0'].value,result.params['E0'].value,result.params['W'].value,result.params['s'].value))
        probb.append(float(stuff[1]))
    else:
        sigmas.append(0.)
        probb.append(float(stuff[1]))

  SEU_xs=0
  for i in range(len(sigmas)):
    if i <len(probb)-1:
      SEU_xs += (sigmas[i+1]-sigmas[i])*probb[i]*1e+08 # probb is the SEU prob in a sensitive volume of 1 cubic micrometer, the cross section is on 1 square micrometer, so 10-8 cm-2

  return SEU_xs     
      
# Compute probability for a given cross section and E0
# Basically like if the Weibull was an Heaviside with a plateau at sigma for
# a threshold of E0
# This useful with no observation. Here your define an upper limit for the Xsec, and use it as sigma.

def SEU_xsec_lim(thresh,sigma):
    sigmas = []
    probb = []

    for row in lineprob:
        stuff=row.split('\n')[0].split(',')
        if float(stuff[0])>thresh:
            sigmas.append(sigma)
        else:
            sigmas.append(0)
        probb.append(float(stuff[1]))
            
    SEU_xs=0
    for i in range(len(sigmas)):
        #print(sigmas[i], probb[i])
        if i <len(probb)-1:
            SEU_xs += (sigmas[i+1]-sigmas[i])*probb[i]*1e+08
                            
    return SEU_xs

'''
The standard macro starts here
'''

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:m:t:",["file=","mode=","type="])
except getopt.GetoptError:
    outlog.write('UpsetRate.py -f <FileName> -m <Mode> -t <Type>')
    sys.exit(2)

file=''
mode=''
typ=''
verbose=0

for opt, arg in opts:
    if opt == '-h':
        outlog.write('UpsetRate.py -f <FileName> -m <Mode> -t <Type>')
        sys.exit()
    elif opt in ("-f", "--file"):
        print(arg)
        file=arg
    elif opt in ("-m", "--mode"):
        mode=arg
        print(arg)
    elif opt in ("-t", "--type"):
        typ=arg
        print(arg)


if typ == 'L1':
  ifL1  = True
  rowID    = 11
elif typ == 'TRG':
  ifL1  = False
  rowID    = 7
else:
  print("The type of data must be L1 or TRG")
  sys.exit()

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

print(file)
fsplit=file.split('/')
print(fsplit)
fname=fsplit[len(fsplit)-1]
print(fname)

datf=open(fname,'r')
probf=open('SEUrate.txt','r')

linedat=datf.readlines()
lineprob=probf.readlines()

data=[]
prob=[]

# Retrieve the data

'''
row id meaning (in file results.txt)

row[0] : data type (MPA of CBC)
row[1] : L1 rate (in kHz)
row[2] : run number during the SEU campaign
row[3] : ion name
row[4] : LET of the ion
row[5] : ion flux, in particles/s/cm2
row[6] : run duration (in s)
row[7] : number of stub simple errors
row[11]: number of L1 simple errors
'''

for line in linedat:
    if line.find(mode)!=-1:
        data.append(line.split(','))

for line in lineprob:
    prob.append(line.split(','))

titles= []
names = []

print("\n")
print("######################################################")
print("################### Reading",data,"data ##############")
print("######################################################")
print("\n")

countHO=0
countLO=0
erry_up=0
erry_down=0
conflevel=0.95
alpha=1-conflevel

xsecdat=[]
xsecdat_up=[]
xsecdat_down=[]
letdat=[]
edepdat=[]

# First read the data corresponding to selected data for selected mode

for row in data:
    if verbose==1:
        print(row)
    nevt=float(row[rowID])
    
    if nevt!=0  and row[1]!='10' and float(row[4])<200:
        
        # First compute the errors
        #
        # nevt is small so error is usually poissonian
        # We want to measure the 95% error bars which means:
        # nevt/+up/-down with:
        #
        # P(obs>up)<2.5%
        # P(obs<down)<2.5%
        #
        # So 95% of the observed values are between up and down.
        # This is given by the function gammaincinv
        #
        # down  = gammaincinv(nevt,0.025)
        # up    = gammaincinv(nevt,0.975)

        down=sp.special.gammaincinv(nevt,(1-conflevel)/2)
        up=sp.special.gammaincinv(nevt,(1+conflevel)/2)
        
        if verbose==1:
           print("down limit = ",down)
           print("up limit   = ",up)
           print(nevt,"events were observed, the 95% CL range of observations is:")
           print(down," < nevt < ",up)
        
        erry_down=(nevt-down)/nevt
        erry_up=(up-nevt)/nevt
        
        # Then compute the cross section
        fluence=float(row[5])*float(row[6])
        xsec = nevt/fluence
        
        # One adds 5% error on fluence measurement
        
        erryu = xsec*math.sqrt(erry_up*erry_up+0.0025)
        erryd = xsec*math.sqrt(erry_down*erry_down+0.0025)
        erry=max(erryu,erryd)

        #print("Cross-section at LET= "+float(row[4])+": ",xsec,"+",erryu,"/-",erryd)
        
        xsecdat.append(xsec)
        xsecdat_up.append(xsec+erryu)
        xsecdat_down.append(xsec-erryd)
        letdat.append(float(row[4]))
        edepdat.append(float(row[4])*0.232)


asymmetric_error = [np.subtract(xsecdat,xsecdat_down), np.subtract(xsecdat_up,xsecdat)]

plt.errorbar(letdat, xsecdat, yerr=asymmetric_error, fmt='ko', label='measurements')
'''
plt.plot(letdat, xsecdat, 'k-', label='measurements')
plt.plot(letdat, xsecdat_up,'-r', label='upper 95% limit')
plt.plot(letdat, xsecdat_down,'-r', label='lower 95% limit')
plt.fill_between(letdat, xsecdat_up, xsecdat_down, color="k", alpha=0.15)
plt.legend()
'''
plt.xlabel("LET (in MeV/$cm^2$/mg)")
plt.ylabel("$\sigma$ (in $cm^2$)")
#plt.xscale("log")
plt.yscale("log")
plt.grid(True, which="both", ls="-")

# create a set of Parameters
params = Parameters()
parnames = ["sigma0","E0","W", "s"]
params.add(parnames[0], value=1e-5, min=1e-7)
params.add(parnames[1], value=0.1,min=0,max=0.5)
params.add(parnames[2], value=2.0, min=0, max=50)
params.add(parnames[3], value=.5, min=0, max=1)

# do fit, here with the default leastsq algorithm
minner = Minimizer(fcn2min,  params, fcn_args=(letdat, xsecdat))
result = minner.minimize()
final = xsecdat + result.residual
# write error report
report_fit(result)

paramsE = Parameters()
paramsE.add(parnames[0], value=1e-5, min=1e-7)
paramsE.add(parnames[1], value=0.02,min=0,max=0.1)
paramsE.add(parnames[2], value=2.0, min=0, max=50)
paramsE.add(parnames[3], value=.5, min=0, max=1)
minner = Minimizer(fcn2min, paramsE, fcn_args=(edepdat, xsecdat))
resultE = minner.minimize()
report_fit(resultE)

minner = Minimizer(fcn2min, paramsE, fcn_args=(edepdat, xsecdat_down))
resultEdown = minner.minimize()
minner = Minimizer(fcn2min, paramsE, fcn_args=(edepdat, xsecdat_up))
resultEup = minner.minimize()

plt.plot(letdat, final,'r--',label='fit')
plt.legend()

xsec_conv=SEU_xsec(resultE)
xsec_convd=xsec_conv-SEU_xsec(resultEdown)
xsec_convu=SEU_xsec(resultEup)-xsec_conv


print("SEU xs for",mode,"case: ", "{:0.2e}".format(xsec_conv))
print("-> Error rate per on ",typ," path:", "{:0.2e}".format(xsec_conv*hrate),"+","{:0.2e}".format(xsec_convu*hrate),"-","{:0.2e}".format(xsec_convd*hrate),"error/s")
print("-> For all modules: ", "{:0.3}".format(xsec_conv*hrate*2*nmods),"+","{:0.3}".format(xsec_convu*hrate*2*nmods),"-","{:0.3}".format(xsec_convd*hrate*2*nmods),"error(s)/s")

plt.savefig("X_section_SEU"+mode+"_"+typ+".png")

