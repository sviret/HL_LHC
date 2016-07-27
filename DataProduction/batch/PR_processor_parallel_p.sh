#!/bin/bash
#
# This script is the main interface for pattern recognition on
# CMSSW files.
#
# It is called by AMPR_parallel_XRD.csh

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# You're not supposed to touch anything here
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if [ ${1} = "PROC" ]; then

#    echo "$PACKDIR/PR_processor_parallel.sh PROC $INDIR_XROOT/$l $OUTD -1 $OUTDIR_GRID $RELEASEDIR $GTAG $OUTDIRSTO" >> run_FIT_${nsj}_${MATTER}.sh

    INPUT=${2}                # The input xrootd file name and address
    OUTPUT=${3}               # Output file name 
    NEVT=${4}                 # #evts/file
    OUTDIR=${5}               # The first event to process in the input file
    CMSSW_PROJECT_SRC=${6}    # The CMSSW project release dir
    GT=${7}                   # The global tag
    INTMP=${8}               # 

    INFILE=`basename $INPUT`
    echo $INPUT,$INFILE
    echo $PWD

    #
    # Setting up environment variables
    #   

    mkdir $INTMP
    cd $INTMP
    TOP=$PWD

    mkdir RECOVERY

    #
    # And we tweak the python generation script according to our needs
    #  

    cd $CMSSW_PROJECT_SRC
    export SCRAM_ARCH=slc6_amd64_gcc472
    eval `scramv1 runtime -sh`   


    cd $TOP
#    cp $CMSSW_PROJECT_SRC/src/DataProduction/test/base/SLHC_EXTR_TILTED_BASE.py  PROC_dummy_${OUTPUT}.py 
    cp $CMSSW_PROJECT_SRC/src/DataProduction/test/base/SLHC_EXTR_FLAT_BASE.py  PROC_dummy_${OUTPUT}.py 

    # Finally the script is modified according to the requests

#    sed "s/NEVTS/$NEVT/"                                   -i PROC_dummy_${OUTPUT}.py
#    sed "s#INPUTFILENAME#file:$INFILE#"                    -i PROC_dummy_${OUTPUT}.py
    sed "s#INPUTFILENAME#$INPUT#"                          -i PROC_dummy_${OUTPUT}.py
    sed "s#OUTFILENAME#$OUTPUT#"                           -i PROC_dummy_${OUTPUT}.py
    sed "s/MYGLOBALTAG/$GT/"                               -i PROC_dummy_${OUTPUT}.py

    cmsRun PROC_dummy_${OUTPUT}.py 

    # Recover the data
    #  

    mv $TOP/$OUTPUT $TOP/RECOVERY/$OUTPUT 

    rm *_dummy*${OUTPUT}.py 

fi

