#!/bin/bash

###########################################
#
# Main script for data generation with tracker only
#
# The jobs themselves are launched by generator_SLHC.sh
#
# Usage:
# source launch_SLHC.sh p1 p2 p3 p4 p5 p6 p7 ...
# with:
# p1 : global tag for MC production (eg auto:upgradePLS3)
# p2 : the number of jobs
# p3 : the number of events per run 
# p4 : MU/PI/ELE/TT/JETS/MUBANK/ELEBANK
# p5 : PTMIN    
# p6 : PTMAX
# p7 : PHIMIN
# p8 : PHIMAX
# p9 : ETAMIN
# p10: ETAMAX
# p11: TOWER ID: the trigger tower ID (not implemented anymore, put -1)
# p12: THRESH: the pt threshold for the stub maker, in GeV (2/3/4, -1 for baseline, 140/200/300 for SW140/200/300 tunings, or 0 for no thresh (FE hardware limits))
# p13: SUFFIX: a specific nickname for this production 
# p14: GEOMTYPE: TILT or FLAT 
# p15: PU: PU scale (0,140,.... actually tested up to 400, but can certainly go a step further)
# p16: BATCH or nothing: launch lxbatch jobs or not
#
# For more details, and examples, have a look at:
# 
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto920 (STEP II)
# https://www.dropbox.com/s/nxajjbwgmnyc0a5/Getting_tilted_stubs_10.pdf?dl=0
# 
# Author: S.Viret (s.viret_at_ipnl.in2p3.fr)
# Date  : 21/12/2017
#
# Script tested with release CMSSW_10_0_0_pre1
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

# The queue over which you want to send the job
BQUEUE=1nd

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

voms-proxy-init --voms cms --valid 100:00 -out $HOME/.globus/gridproxy.cert
export X509_USER_PROXY=${HOME}'/.globus/gridproxy.cert'

GTAG=${1}           # Global tag
N_RUN=${2}          # Number of samples 
EVTS_PER_RUN=${3}   # Number of events per sample
MATTER=${4}         # Type of event
GTYPE=${14}         # Type of geometry
PTYPE=13            # Default particle type (MUON)
GEOM=False          # Default geometry is tilted
NPU=${15}           # 
i=0

# The geomtry type is set depending on the GTYPE type
case "$GTYPE" in
"TILT")
GEOM=False;
;;
"FLAT")
GEOM=True;
;;
*)
GTYPE=FLAT;
echo "You have not provided a correct geom name, we will use FLAT by default";
;;
esac

# The particle type is set depending on the MATTER type

# First, single particle guns

case "$MATTER" in
"MU")
PTYPE=13;
;;
"PI")
PTYPE=211;
;;
"ELE")
PTYPE=11;
;;
"JETS")
PTYPE=1;
;;
"TT")
PTYPE=2;
;;
"MUBANK")
PTYPE=1013;
;;
"ELEBANK")
PTYPE=1011;
;;
"SWTUNING")
PTYPE=666;
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

DIRNAME=${MATTER}_PU${NPU}_${GTYPE}_${13}

# Then we setup some environment variables

cd  ..
PACKDIR=$PWD           # This is where the package is installed 
cd  ../..
RELEASEDIR=$PWD        # This is where the release is installed

cd $PACKDIR/batch

# Finally we create the batch scripts and launch the jobs

echo 'Creating directory '$STORAGEDIR/$DIRNAME' for type '$PTYPE

OUTPUTDIR=$STORAGEDIR2/$DIRNAME
OUTPUTDIR2=$STORAGEDIR/$DIRNAME

gfal-mkdir $OUTPUTDIR

# Special case, rerun the stub production on relval sample

if [ $MATTER == "SWTUNING" ]; then

    dasgoclient --query="dataset=/RelVal"${13}"/${CMSSW_VERSION}_*_upgrade2023_realistic*2023D21*"$NPU"*/GEN-SIM-DIGI-RAW" --limit=1 | grep GEN-SIM > temp

    samplelist=`cat temp`


    for f in $samplelist
    do

	type=`echo $f | cut -d'/' -f2`
	rel=`echo $f | cut -d'/' -f3 | cut -d'-' -f1`
	name=`echo $f | cut -d'/' -f3`

	echo $type,$rel,$name

	echo " #!/bin/bash " > files.txt
	dasgoclient --query="file dataset="$f --limit=50 | grep store >> files.txt

	filelist=`cat files.txt`

	for f in $filelist
	do

	    echo "#\!/bin/bash" > gen_job_${DIRNAME}_${i}.sh
	    echo "source $PACKDIR/batch/generator_SLHC.sh $EVTS_PER_RUN $PTYPE $MATTER $GTAG $i $RELEASEDIR $PACKDIR $OUTPUTDIR ${PTMIN} ${PTMAX} ${PHIMIN} ${PHIMAX} ${ETAMIN} ${ETAMAX} $NPU $THRESH $TID $GEOM $f" >> gen_job_${DIRNAME}_${i}.sh
	    chmod 755 gen_job_${DIRNAME}_${i}.sh

	    if [ ${16} == "BATCH" ]; then	
		bsub -q $BQUEUE -e /dev/null -o /tmp/${LOGNAME}_out.txt gen_job_${DIRNAME}_${i}.sh
	    fi

	    i=$(( $i + 1))

	done
    done   
fi  

# Normal case

while [ $i -lt $N_RUN ]
do
    rm *.log

    if [ $MATTER == $MATTER ]; then

	deal=`gfal-ls $OUTPUTDIR/SLHC_extr_PU${NPU}_${MATTER}_${i}.root | wc -l`

	if [ $deal != "0" ]; then

	    i=$(( $i + 1))
	    continue
	fi
    fi 

    echo "#\!/bin/bash" > gen_job_${DIRNAME}_${i}.sh
    echo "source $PACKDIR/batch/generator_SLHC.sh $EVTS_PER_RUN $PTYPE $MATTER $GTAG $i $RELEASEDIR $PACKDIR $OUTPUTDIR ${PTMIN} ${PTMAX} ${PHIMIN} ${PHIMAX} ${ETAMIN} ${ETAMAX} $NPU $THRESH $TID $GEOM" >> gen_job_${DIRNAME}_${i}.sh
    chmod 755 gen_job_${DIRNAME}_${i}.sh

    if [ ${16} == "BATCH" ]; then	
	bsub -q $BQUEUE -e /dev/null -o /tmp/${LOGNAME}_out.txt gen_job_${DIRNAME}_${i}.sh
    fi

    i=$(( $i + 1))
done


