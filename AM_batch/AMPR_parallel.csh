#!/bin/csh

###########################################
#
# Main script for parallel CMSSW AM-based pattern 
# recognition in batch on a multi-core machine
# 
# !!!Working on EDM files!!!
#
# !!!Requires the installation of GNU parallel executable on your machine!!!
#    !!!!!!! http://www.gnu.org/software/parallel/ !!!!!!!
#
# The jobs themselves are launched by PR_processor_parallel.sh
#
# source AMPR.csh p1 p2 p3 p4 p5 p6
# with:
# p1 : The subdirectory containing the data file you want to analyze (best is to copy them beforehand on the machine scratch area)
# p2 : The directory where you will retrieve the bank files, the pattern reco will
#      run over all the pbk files contained in this directory
# p3 : How many events per input data file? 
# p4 : How many events per job (should be below p3...)?
# p5 : The global tag name
# p6 : How many cores you want to use in parallel
#
# For more details, and examples, have a look at:
# 
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP V.2)
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 28/04/2014
#
# Script tested with release CMSSW_6_2_0_SLHC10
#
###########################################


# Here we retrieve the main parameters for the job 

set MATTER  = ${1}   # Directory where the input root files are
set BANKDIR = ${2}   # Directory where the bank (.pbk) files are
set NTOT    = ${3}   # How many events per data file?
set NPFILE  = ${4}   # How many events per job?
set GTAG    = ${5}   # Global tag
set NCORES  = ${6}   # #cores
set NFILES  = ${7}   # #files per job

###################################
#
# The list of parameters you can modify is here
#
###################################

# You have to adapt this to your situation

# The scratch directory where you put the input and temporary files
# !!! Ensure that you have enough scratch space available !!!

set BASE = /scratch/viret

set INDIR       = $BASE/$MATTER
mkdir $BASE/TMP
set OUTDIRTMP   = $BASE/TMP/$MATTER

# The SE directory containing the output EDM file with the PR output
set OUTDIR      = /dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/PR/$MATTER

setenv LFC_HOST 'lyogrid06.in2p3.fr'
set INDIR_GRID   = srm://$LFC_HOST/$INDIR
set INDIR_XROOT  = root://$LFC_HOST/$INDIR
set OUTDIR_GRID  = srm://$LFC_HOST/$OUTDIR
set OUTDIR_XROOT = root://$LFC_HOST/$OUTDIR

# The parallel command
set parallel = /gridgroup/cms/brochet/.local/bin/parallel

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
setenv X509_USER_PROXY ${HOME}'/.globus/gridproxy.cert'


# Then we setup some environment variables

cd  ..
set PACKDIR      = $PWD           # This is where the package is installed 
cd  ../..
set RELEASEDIR   = $CMSSW_BASE    # This is where the release is installed

cd $PACKDIR/batch


echo 'The data will be read from directory: '$INDIR
echo 'The final pattern reco output files will be written in: '$OUTDIR

lfc-mkdir $OUTDIR
mkdir $OUTDIRTMP
mkdir ${INDIR}/TMP

# We loop over the data directory in order to find all the files to process

@ ninput = 0	 
@ nsj    = 0
@ npsc   = $NFILES

echo "#\!/bin/bash" > global_stuff.sh

foreach ll (`\ls $INDIR | grep EDM`) 
   
    set l = `basename $ll`

    @ i = 0
    @ j = $NPFILE

    set val = `expr $ninput % $npsc`

    if ($val == 0) then
	@ nsj += 1

	echo "source $PACKDIR/batch/run_${nsj}.sh"  >> global_stuff.sh

	echo "#\!/bin/bash" > run_${nsj}.sh
	echo "${parallel} -j ${NCORES} < $PACKDIR/batch/run_PR_${nsj}.sh" >> run_${nsj}.sh
	echo "${parallel} -j ${NCORES} < $PACKDIR/batch/run_MERGE_${nsj}.sh" >> run_${nsj}.sh
	echo "${parallel} -j ${NCORES} < $PACKDIR/batch/run_FMERGE_${nsj}.sh" >> run_${nsj}.sh
	echo "${parallel} -j ${NCORES} < $PACKDIR/batch/run_FIT_${nsj}.sh" >> run_${nsj}.sh
	echo "${parallel} -j ${NCORES} < $PACKDIR/batch/run_RM_${nsj}.sh" >> run_${nsj}.sh

	echo "#\!/bin/bash" > run_PR_${nsj}.sh
	echo "#\!/bin/bash" > run_MERGE_${nsj}.sh
	echo "#\!/bin/bash" > run_FMERGE_${nsj}.sh
	echo "#\!/bin/bash" > run_FIT_${nsj}.sh
	echo "#\!/bin/bash" > run_RM_${nsj}.sh

	chmod 755 run_${nsj}.sh
	chmod 755 run_PR_${nsj}.sh
	chmod 755 run_MERGE_${nsj}.sh
	chmod 755 run_FMERGE_${nsj}.sh
	chmod 755 run_FIT_${nsj}.sh
	chmod 755 run_RM_${nsj}.sh

    endif

    @ ninput += 1

    echo 'Working with file '$l

    # First look if the file has been processed

    set OUTM  = `echo $l | cut -d. -f1`

    set OUTF  = ${OUTM}"_with_AMPR.root"
    set OUTE  = ${OUTM}"_with_FIT.root"
    set OUTD  = ${OUTM}"_extr.root"

    @ processed = 0
    @ section   = 0

    while ($i < $NTOT)
		
	@ sec = 0
        @ secdone = 0
        @ section += 1

	#
	# First step, we loop over the banks and run the 
	# AM PR on the given data sample
	#

        foreach k (`\ls $BANKDIR | grep _sec`) 

	    # By default, for CMSSW, we loop over all available bank in the directory provided

	    set SECNUM = `echo $k | cut -d_ -f5 | sed s/sec//` 
	    set OUTS1 = `echo $l | cut -d. -f1`_`echo $k | cut -d. -f1`_${i}_${j}

	    echo "$PACKDIR/batch/PR_processor_parallel.sh PR ${INDIR}/$l $BANKDIR/$k $OUTS1.root  ${i} $NPFILE $OUTDIRTMP $RELEASEDIR $sec $GTAG $SECNUM ${INDIR}/TMP" >> run_PR_${nsj}.sh
	    @ sec += 1

	end # End of bank loop

	#
	# Second step, for this given part of the file, all the  
	# banks output are available. We then launch the merging 
	# procedure
	#

	echo "$PACKDIR/batch/PR_processor_parallel.sh  MERGE ${i}_${j}.root $OUTDIRTMP $OUTDIRTMP ${OUTM}_ $RELEASEDIR $GTAG ${OUTM}_${i}_${j}" >> run_MERGE_${nsj}.sh

	@ i += $NPFILE
	@ j += $NPFILE

    end # End of loop over one input file
	
    #
    # Third step, all the merged files for the given input
    # file have been processed. Then launch the final merging 
    # 

    echo "$PACKDIR/batch/PR_processor_parallel.sh  FINAL MERGED_${OUTM}_ $OUTDIRTMP $OUTDIRTMP $OUTF $RELEASEDIR ${OUTM}"  >> run_FMERGE_${nsj}.sh
    echo "$PACKDIR/batch/PR_processor_parallel.sh  FIT $OUTDIRTMP/${OUTF} $OUTE $OUTD $NTOT $OUTDIR_GRID $RELEASEDIR $GTAG ${OUTM}" >> run_FIT_${nsj}.sh
    echo "rm $OUTDIRTMP/*${OUTM}_*" >> run_RM_${nsj}.sh

end

chmod 755 global_stuff.sh

