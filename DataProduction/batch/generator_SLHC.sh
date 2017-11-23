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
NPU=${15}               # Average PU events number
THRESH=${16}            # pt threshold for stub windows
TID=${17}               # trigger tower ID in which the PGun particles are sent
GEOM=${18}              # geometry type: FLAT (True) or TILTED (False)
FILE=${19}              # input file name (for SW tuning only)
PU=0 
BANK=0
SW=0

#
# Setting up environment variables
#

TOP=$PWD

cd $CMSSW_PROJECT_SRC
export SCRAM_ARCH=slc6_amd64_gcc530
eval `scramv1 runtime -sh`   
voms-proxy-info


# We take care of random numbers, so that they are different for all runs

MAXCOUNT=$NRUN
START=$(($NRUN * $EVTS)) # To keep track of the right event number (for parallel production)

echo $START

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
SEED5=$RANDOM
SEED6=$RANDOM
SEED7=$RANDOM
SEED8=$RANDOM

echo $PTYPE

#
# And we tweak the python generation script according to our needs
#

cd $TOP

# Download the main config script (Default is particle gun)
cp $PACK_DIR/test/base/SLHC_PU_BASE.py BH_dummy.py 

# Script for Jets production
if [ "$PTYPE" -eq 1 ]; then
    cp $PACK_DIR/test/base/SLHC_JET_BASE.py BH_dummy.py 
fi

# Script for TTbar production
if [ "$PTYPE" -eq 2 ]; then
    cp $PACK_DIR/test/base/SLHC_TT_BASE.py BH_dummy.py 
fi

cp $PACK_DIR/test/base/SLHC_MBIAS_BASE.py BH_dummyMinBias.py 
cp $PACK_DIR/test/base/SLHC_EXTR_BASE.py BH_dummy2.py 


# Scripts for AM bank sample production

if [ "$PTYPE" -eq 1013 ]; then
    PTYPE=13
    BANK=1
    cp $PACK_DIR/test/base/SLHC_BANKSIM_BASE.py BH_dummy.py 
    cp $PACK_DIR/test/base/SLHC_BANK_BASE.py BH_dummy2.py 
fi

if [ "$PTYPE" -eq 1011 ]; then
    PTYPE=11
    BANK=1
    cp $PACK_DIR/test/base/SLHC_BANKSIM_BASE.py BH_dummy.py 
    cp $PACK_DIR/test/base/SLHC_BANK_BASE.py BH_dummy2.py 
fi

if [ "$PTYPE" -eq 666 ]; then
    SW=1
    cp $PACK_DIR/test/base/SLHC_StubsExtr_BASE.py BH_dummy.py 
fi


# Choice of the stub windows

if [ "$THRESH" -eq 0 ]; then
    echo 'No stub windows!'
    cp $PACK_DIR/test/base/NoTune.txt tune 
elif [ "$THRESH" -eq -1 ]; then
    echo 'Using the default stub windows'
    cp $PACK_DIR/test/base/BaseTune.txt tune 
elif [ "$THRESH" -eq 140 ]; then
    echo 'Using the SW140 stub windows'
    cp $PACK_DIR/test/base/BaseTune_SW140.txt tune 
elif [ "$THRESH" -eq 200 ]; then
    echo 'Using the SW200 stub windows'
    cp $PACK_DIR/test/base/BaseTune_SW200.txt tune 
else
    cp $PACK_DIR/test/base/${THRESH}GevTune.txt tune 
fi

SWTUN=`cat tune`

echo $SWTUN
echo $FILE

# Finally the script is modified according to the requests

sed "s/TILTEDORFLAT/$GEOM/"                            -i BH_dummy.py
sed "s/NEVTS/$EVTS/"                                   -i BH_dummy.py
sed "s/RUN/$START/"                                    -i BH_dummy.py
sed "s/PTYPE/$PTYPE/"                                  -i BH_dummy.py
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
sed "s#PUFILEA#$TOP/MBiasSample.root#"                 -i BH_dummy.py

sed "s/TILTEDORFLAT/$GEOM/"                            -i BH_dummy2.py
sed "s/MYGLOBALTAG/$GTAG/"                             -i BH_dummy2.py

perl -i.bak -pe 's/SWTUNING/'"${SWTUN}"'/g' BH_dummy.py

if [ "$BANK" -eq 1 ]; then
    sed "s/NEVTS/$EVTS/"                                   -i BH_dummy2.py
    sed "s/MYGLOBALTAG/$GTAG/"                             -i BH_dummy2.py
fi


EVTS=$(( 5 * $NPU + 1 )) 
sed "s/TILTEDORFLAT/$GEOM/"                            -i BH_dummyMinBias.py
sed "s/NEVTS/$EVTS/"                                   -i BH_dummyMinBias.py
sed "s/NSEEDA/$SEED5/g"                                -i BH_dummyMinBias.py
sed "s/NSEEDB/$SEED6/g"                                -i BH_dummyMinBias.py
sed "s/NSEEDC/$SEED7/g"                                -i BH_dummyMinBias.py
sed "s/MYGLOBALTAG/$GTAG/"                             -i BH_dummyMinBias.py


# Set output filenames
#

DATA_NAME=SLHC_extr_$tag.root

# Launch the whole thing
#

if [ "$BANK" -eq 0 ]; then
    cmsRun BH_dummyMinBias.py
fi

cmsRun BH_dummy.py

if [ "$SW" -eq 0 ]; then
    cmsRun BH_dummy2.py
fi

if [ "$BANK" -eq 1 ]; then
    rm PGun_example.root
    mv extracted_skimmed.root extracted.root
fi

# Recover the data
#

ls -l
gfal-copy file://$TOP/extracted.root $STOR_DIR/$DATA_NAME

if [ "$SW" -eq 0 ]; then
    gfal-copy file://$TOP/PGun_example.root $STOR_DIR/EDM_$DATA_NAME
fi

