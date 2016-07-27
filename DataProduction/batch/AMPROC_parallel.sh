#!/bin/bash

###########################################
#
# Main script for parallel CMSSW fitting stage 
# in batch on a multi-core machine
# 
# !!!Working on EDM files!!!
#
# !!!Requires the installation of GNU parallel executable on your machine!!!
#    !!!!!!! http://www.gnu.org/software/parallel/ !!!!!!!
#
# If you cannot use parallel, set p4 to 1 in the launch command
#
# The jobs themselves are launched by PR_processor_parallel_p.sh
#
# source AMFIT_parallel.sh p1 p2 p3 p4 p5
# with:
# p1 : The directory containing the data file you want to analyze (best is to copy them beforehand on the machine scratch area)
# p2 : Name of the SE subdirectory where you will store the data
# p3 : The global tag name
# p4 : How many cores you want to use in parallel (if one then parallel is not used)
# p5 : How many events per job to process
#
# For more details, and examples, have a look at:
# 
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP V.2)
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 28/04/2014
#
# Script tested with release CMSSW_6_2_0_SLHC14
#
###########################################


# Here we retrieve the main parameters for the job 

INPUTDIR=${1}  # Directory where the input root files are
OUTPUTDIR=${2} # Name of the directory in the SE
GTAG=${3}      # Global tag
NCORES=${4}    # #cores
NFILES=${5}    # #files per job

###################################
#
# The list of parameters you can modify is here
#
###################################


export LFC_HOST=lyogrid06.in2p3.fr

# I/O Directory info on the grid
BAREPATH=/dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/GEN
SRMPATH=srm://$LFC_HOST/$BAREPATH
XRDPATH=root://$LFC_HOST/$BAREPATH

# The parallel command

parallel=/gridsoft/ipnls/parallel/bin/parallel

###########################################################
###########################################################
# You are not supposed to touch the rest of the script !!!!
###########################################################
###########################################################


# The SE directory containing the input EDM files

INDIR_GRID=$SRMPATH/$INPUTDIR
INDIR_XROOT=$XRDPATH/$INPUTDIR

# The SE directory containing the output EDM file with the PR output
OUTDIR_GRID=$SRMPATH/$OUTPUTDIR
OUTDIR_XROOT=$XRDPATH/$OUTPUTDIR

MATTER=`basename $OUTPUTDIR`

OUTDIRTMP=/data/viret/AMProc/TMPDAT_$MATTER
OUTDIRSTO=/data/viret/AMProc/TMP_$MATTER

echo 'The data will be read from directory: '$INDIR_GRID
echo 'The final pattern reco output files will be written in: '$OUTDIR_GRID

# Following lines suppose that you have a certificate installed on lxplus. To do that follow the instructions given here:
#
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideLcgAccess
#

voms-proxy-init --voms cms --valid 100:00 -out $HOME/.globus/gridproxy.cert
export X509_USER_PROXY=${HOME}/.globus/gridproxy.cert

lfc-mkdir -p $BAREPATH/$OUTPUTDIR

# Then we setup some environment variables

PACKDIR=$PWD            # This is where the scripts are created 
RELEASEDIR=$CMSSW_BASE  # This is where the release is installed

mkdir $OUTDIRTMP
mkdir $OUTDIRSTO

# We loop over the data directory in order to find all the files to process

ninput=0	 
nsj=0
npsc=$NFILES

echo "#\!/bin/bash" > global_stuff_${MATTER}.sh

for ll in `lcg-ls $INDIR_GRID | grep EDM` 
do   
    l=`basename $ll`

    val=`expr $ninput % $npsc`

    if [ $val = 0 ]; then

	nsj=$(( $nsj + 1))

	echo "source $PACKDIR/run_${nsj}_${MATTER}.sh"  >> global_stuff_${MATTER}.sh

	echo "#\!/bin/bash" > run_${nsj}_${MATTER}.sh

	if [ $NCORES = 1 ]; then
	    echo "$PACKDIR/run_EXTR_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	    echo "$PACKDIR/run_RM_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	else
	    echo "${parallel} -j ${NCORES} < $PACKDIR/run_EXTR_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	    echo "$PACKDIR/run_RM_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
	fi

	echo "#\!/bin/bash" > run_EXTR_${nsj}_${MATTER}.sh
	echo "#\!/bin/bash" > run_RM_${nsj}_${MATTER}.sh

	chmod 755 run_EXTR_${nsj}_${MATTER}.sh
	chmod 755 run_RM_${nsj}_${MATTER}.sh

    fi

    ninput=$(( $ninput + 1))

#    echo 'Working with file '$l

    # First look if the file has been processed

    OUTB=`echo $l | cut -d. -f1`

    OUTD=${OUTB}"_extracted.root"

    FIT_FILE=${OUTDIRTMP}/RECOVERY/$OUTD

    deal=`ls $FIT_FILE | wc -l`

    if [ $deal != "0" ]; then
	
	echo "File "$l" has already processed, skip it..."
	continue;
    fi
	
    #
    # Last step, all the merged files for the given input
    # file have been processed. Then launch the final merging 
    # 

    echo "$PACKDIR/PR_processor_parallel_p.sh PROC $INDIR_XROOT/$l $OUTD -1 $OUTDIR_GRID $RELEASEDIR $GTAG $OUTDIRTMP" >> run_EXTR_${nsj}_${MATTER}.sh
    echo "cd $OUTDIRTMP" >> run_RM_${nsj}_${MATTER}.sh
    echo "rm *${OUTM}_*root" >> run_RM_${nsj}_${MATTER}.sh
done

chmod 755 global_stuff_${MATTER}.sh

