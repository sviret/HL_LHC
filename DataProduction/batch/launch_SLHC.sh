#!/bin/bash


###########################################
#
# Main script for data generation
#
# The jobs themselves are launched by generator_SLHC.sh
#
# Usage:
# source launch_SLHC.sh p1 p2 p3 p4 p5 p6 p7 ...
# with:
# p1 : global tag for MC production (eg POSTLS161_V15::All)
# p2 : the number of runs
# p3 : the number of events per run 
# p4 : MU/ANTIMU/PIP/PIM/....
# p5 : PTMIN 
# p6 : PTMAX
# p7 : PHIMIN
# p8 : PHIMAX
# p9 : ETAMIN
# p10: ETAMAX
# p11: TOWER ID: the trigger tower ID (put -1 for default)
# p12: THRESH: the pt threshold for the stub maker (2/3/4 of -1 for no thresh)
# p13: SUFFIX: a specific nickname for this production 
# p14: BATCH or nothing: launch lxbatch jobs or not
#
# For more details, and examples, have a look at:
# 
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP II)
# 
# Author: S.Viret (viret@in2p3.fr)
# Date  : 12/04/2013
#
# Script tested with release CMSSW_6_2_0_SLHC14
#
###########################################



###################################
#
# The list of parameters you can modify is here
#
###################################

# You have to adapt this to the storage element you plan to use

# Info concerning the grid directory where the data will be stored
export LFC_HOST=lyogrid06.in2p3.fr
STORAGEDIR=/dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/GEN 
STORAGEPU=yourdir
# The average number of PU event (for PU simulation)
NPU=140

# The queue over which you want to send the job
BQUEUE=8nh

# The directory GRID-friendly address
STORAGEDIR2=srm://$LFC_HOST/$STORAGEDIR                            

###########################################################
###########################################################
# You are not supposed to touch the rest of the script !!!!
###########################################################
###########################################################

# Following lines suppose that you have a certificate installed on lxplus. To do that follow the instructions given here:
#
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideLcgAccess
#

source /afs/cern.ch/project/gd/LCG-share/current_3.2/etc/profile.d/grid-env.sh
voms-proxy-init --voms cms --valid 100:00 -out $HOME/.globus/gridproxy.cert
export X509_USER_PROXY=${HOME}'/.globus/gridproxy.cert'

GTAG=${1}           # Global tag
N_RUN=${2}          # Number of samples 
EVTS_PER_RUN=${3}   # Number of events per sample
MATTER=${4}         # Type of event
PTYPE=-13           # Default particle type (MUON)

i=0

# The particle type is set depending on the MATTER type

# First, single particle guns

case "$MATTER" in
"MU")
PTYPE=-13;
;;
"ANTIMU")
PTYPE=13;
;;
"PIP")
PTYPE=211;
;;
"PIM")
PTYPE=-211;
;;
"ELE")
PTYPE=-11;
;;
"POS")
PTYPE=11;
;;
"JETS")
PTYPE=1;
;;
"MINBIAS")
PTYPE=666;
;;
"PILEUP")
PTYPE=777;
;;
"PILEUP_ND")
PTYPE=778;
;;
"PILEUP4T")
PTYPE=888;
;;
"PILEUPREDO")
PTYPE=7777;
;;
"MUBANK")
PTYPE=1013;
;;
*)
PTYPE=13; 
;;
esac

# Parameters for particle gun

PTMIN=${5}
PTMAX=${6}
PHIMIN=${7}
PHIMAX=${8}
ETAMIN=${9}
ETAMAX=${10}
TID=${11}
THRESH=${12}

# Then we setup some environment variables

cd  ..
PACKDIR=$PWD           # This is where the package is installed 
cd  ../..
RELEASEDIR=$PWD        # This is where the release is installed

cd $PACKDIR/batch

# Finally we create the batch scripts and launch the jobs

echo 'Creating directory '$STORAGEDIR/${MATTER}_${13}' for type '$PTYPE

OUTPUTDIR=$STORAGEDIR2/${MATTER}_${13}
OUTPUTDIR2=$STORAGEDIR/${MATTER}_${13}

lfc-mkdir $OUTPUTDIR2 

while [ $i -lt $N_RUN ]
do
    rm *.log

    if [ $MATTER == "MUBANK" ]; then

	deal=`lcg-ls $OUTPUTDIR/SLHC_extr_MUBANK_${i}.root | wc -l`

	if [ $deal != "0" ]; then
	
	    run =`ls gen_job_${MATTER}_${12}_${13}_${i}.sh | wc -l`

	    if [ $run != "0" ]; then
		rm gen_job_${MATTER}_${12}_${13}_${i}.sh
	    fi
	    i=$(( $i + 1))
	    continue
	fi

	run=`ls gen_job_${MATTER}_${12}_${13}_${i}.sh | wc -l`

	if [ $run != "0" ]; then
	    i=$(( $i + 1))
	    continue
	fi
    fi 

    echo "#\!/bin/bash" > gen_job_${MATTER}_${12}_${13}_${i}.sh
    echo "source $PACKDIR/batch/generator_SLHC.sh $EVTS_PER_RUN $PTYPE $MATTER $GTAG $i $RELEASEDIR $PACKDIR $OUTPUTDIR ${PTMIN} ${PTMAX} ${PHIMIN} ${PHIMAX} ${ETAMIN} ${ETAMAX} $STORAGEPU $NPU ${THRESH} $TID" >> gen_job_${MATTER}_${12}_${13}_${i}.sh
    chmod 755 gen_job_${MATTER}_${12}_${13}_${i}.sh

    if [[ -z "$14" ]]; then
	if [ ${14} == "BATCH" ]; then	
	    bsub -q $BQUEUE -e /dev/null -o /tmp/${LOGNAME}_out.txt gen_job_${MATTER}_${12}_${13}_${i}.sh
	fi
    fi

    i=$(( $i + 1))
done


