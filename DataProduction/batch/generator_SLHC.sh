#!/bin/bash
#
# This script is the main interface for event production
#
# It is called by launch_SLHC.csh

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# You're not supposed to touch anything here
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#
# First we retrieve the input variables
#

EVTS=${1}               # Number of events to produce
PTYPE=${2}              # Production type
TYPE=${3}               # "
GTAG=${4}               # Global tag
NRUN=${5}               # Run number
CMSSW_PROJECT_SRC=${6}  # CMSSW release dir
PACK_DIR=${7}           # Batch package dir 
STOR_DIR=${8}           # Storage dir
PTMIN=${9}              #
PTMAX=${10}             # 
PHIMIN=${11}            # ==> The particle gun parameters
PHIMAX=${12}            #
ETAMIN=${13}            #
ETAMAX=${14}            #
PUDIR=${15}             # Minbias storage repository
NPU=${16}               # PU events number


#
# Setting up environment variables
#

TOP=$PWD

cd $CMSSW_PROJECT_SRC
export SCRAM_ARCH=slc5_amd64_gcc434
eval `scramv1 runtime -sh`   
voms-proxy-info


# We take care of random numbers, so that they are different for all runs

MAXCOUNT=$NRUN
count=0
tag=${TYPE}_${NRUN}

while [ "$count" -le $MAXCOUNT ]      # Generate 10 ($MAXCOUNT) random integers.
do
  echo $RANDOM
  let "count += 1"  # Increment count.
done

SEED1=$RANDOM
SEED2=$RANDOM
SEED3=$RANDOM
SEED4=$RANDOM


#
# And we tweak the python generation script according to our needs
#

cd $TOP

# Download the main config script
cp $PACK_DIR/test/base/SLHC_PGUN_BASE.py BH_dummy.py 


# Script for MinBias production
if [ "$PTYPE" -eq 666 ]; then
    cp $PACK_DIR/test/base/SLHC_MBIAS_BASE.py BH_dummy.py 
fi

# Scripts for Track trigger production

if [ "$PTYPE" -eq 1013 ]; then
    PTYPE=13
    cp $PACK_DIR/test/base/SLHC_BANK_BASE.py BH_dummy.py 
fi

if [ "$PTYPE" -eq -1013 ]; then
    PTYPE=-13
    cp $PACK_DIR/test/base/SLHC_BANK_BASE.py BH_dummy.py 
fi


# Script for Pileup production

# Pile up requires a bit of tweaking
# Indeed we choose two random MinBias files in which the PU events will be chosen
# They are downloaded on the batch machine
#

FILE1=`echo $(( $RANDOM % 200 ))`
FILE2=`echo $(( $RANDOM % 200 ))`

if [ "$PTYPE" -eq 777 ]; then
    PTYPE=-13
    cp $PACK_DIR/test/base/SLHC_PU_BASE.py BH_dummy.py 
    lcg-cp $PUDIR/SLHC_extr_MINBIAS_$FILE1.root file:/$TOP/SLHC_extr_MINBIAS_$FILE1.root 
    lcg-cp $PUDIR/SLHC_extr_MINBIAS_$FILE2.root file:/$TOP/SLHC_extr_MINBIAS_$FILE2.root 
    tag=${TYPE}_${NPU}_${NRUN}
fi



# Finally the script is modified according to the requests

sed "s/NEVTS/$EVTS/"                                   -i BH_dummy.py
sed "s/PTYPE/$PTYPE/"                                  -i BH_dummy.py
sed "s#INPUTFILENAME#$TOP/SLHC_extr_PU_${NRUN}.root#"  -i BH_dummy.py
sed "s/NSEEDA/$SEED1/g"                                -i BH_dummy.py
sed "s/NSEEDB/$SEED2/g"                                -i BH_dummy.py
sed "s/NSEEDC/$SEED3/g"                                -i BH_dummy.py
sed "s/NSEEDD/$SEED4/g"                                -i BH_dummy.py
sed "s/MYGLOBALTAG/$GTAG/"                             -i BH_dummy.py
sed "s/PTMIN/$PTMIN/"                                  -i BH_dummy.py
sed "s/PTMAX/$PTMAX/"                                  -i BH_dummy.py
sed "s/PHIMIN/$PHIMIN/"                                -i BH_dummy.py
sed "s/PHIMAX/$PHIMAX/"                                -i BH_dummy.py
sed "s/ETAMIN/$ETAMIN/"                                -i BH_dummy.py
sed "s/ETAMAX/$ETAMAX/"                                -i BH_dummy.py
sed "s/NPU/$NPU/"                                      -i BH_dummy.py
sed "s#PUFILEA#$TOP/SLHC_extr_MINBIAS_$FILE1.root#"    -i BH_dummy.py
sed "s#PUFILEB#$TOP/SLHC_extr_MINBIAS_$FILE2.root#"    -i BH_dummy.py



# Set output filenames
#

DATA_NAME=SLHC_extr_$tag.root

# Launch the whole thing
#

cmsRun BH_dummy.py


# Recover the data
#

ls -l
lcg-cp file://$TOP/extracted.root $STOR_DIR/$DATA_NAME

