'''

DataParse_trigger.py

Main data parser for CIC SEU data (TRIGGER part).

Takes as input a txt file containing the info for one run
Provides in output different type of info on TRG paths errors

Input file is the one produced by SEU_test_run_from_cache.py
The bin file used to run the test is also necessary (for the comparison)
It's VERY IMPORTANT to use the exact same bin file here, make sure you have it correct

This file is a first pass, analyzing the most common errors.
More complicated case must be analyzed individually

Author: S.Viret
Date  : 2020

'''


import sys, getopt
import os
import math
import pickle
from datetime import datetime

verbose=0  # Debug mode, put to 1 to get massive outlog.writeouts

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:",["file="])
except getopt.GetoptError:
    outlog.write('DataParse_trigger.py -f <FileName>')
    sys.exit(2)

file=''

for opt, arg in opts:
    if opt == '-h':
        outlog.write('DataParse_trigger.py -f <FileName>')
        sys.exit()
    elif opt in ("-f", "--file"):
        file=arg

#
# STEP 1. Open the file and parse the main info
#

freq=640  # The default output freq is 640Mz

fsplit=file.split('/')
fname=fsplit[len(fsplit)-1]

outputlog='TRG_SEU_ana_'+fname
outlog=open(outputlog,'w')

fp = open(file, 'r')

lines = fp.readlines()

lOutTriggerBitstream=[]
llData=[]
TRGO=[]

errlist=[]
errlist_trg=[]

bitflip_hard=0
bitflip=0
bitflip_trg=0
MPA=1
t2=-1
runtime=0

ncycles=2  # this is the default value, but corrected if necessary

#We loop over all the lines of the output file, and check for data

for line in lines:

    if line.find('elapsed')!=-1:
        runtime=(line.split(':')[2]) # Run lenght (in s)
    
    # Retrieve the reference data bin file injected into the RAM
    # And store the info to perform further comparisons
    #
    # This line is at the start of the file so that we can
    # Retrieve the input frames at the beginning
    #

    if line.find('Standard')!=-1:
        fname=(line.split('file ')[1]).split('.bin')[0]+'.bin'
        if line.find('MPA')==-1:
            MPA=0
        if line.find('_640_')==-1:
            freq=320

        inputbare=fname.split("/")
        filename=inputbare[len(inputbare)-1]
        # filename is the name of the .bin file, in contains the main run
        # characteristics
        inputbare=filename.split(".")
        params=inputbare[0].split("_")
        
        ncycles = int(params[3])
        outlog.write('Opening bin file'+fname+'\n')

        inBitStreamFile = open(fname,'rb')
        input = pickle.load(inBitStreamFile) # Open the bit file w/pickle
        
        lOutTriggerBitstream=input[2]  # The trigger reference output
                                       # This is what we expect from the run
                                       # Into the CIC
        #llData=input[4]                # The output stored on firmware (not used for the moment)

        # TRGO contains the full triggger output over the ncycles generated cycles
        # The 6 first bits contains the first clk cycles (320 or 640) of the 6 output lines, and so on...

        for i in range(len(lOutTriggerBitstream)): # Usually nLines
            word=[]
            for j in range(len(lOutTriggerBitstream[i])):
                word.append(lOutTriggerBitstream[i][j])
            sword="".join(word)
            TRGO.append(sword)
    # End of the input file reading loop

    '''
    # Here we parse the error line provided by the analysis software:
    There are 24 16 bits words containing the following infos:
    
    col 1->6     : CIC trigger lines output
    col 7        : L1 CIC output not realigned with RAM (unused)
    col 8        : Error flags (8 first bits for trigger, 8 following bits for L1).
    col 9->14    : RAM trigger lines (to be compared col 1->6) // The REFERENCE
    col 15       : RAM L1 line // The REFERENCE
    col 16       : CIC L1 output realigned (to be compared with 15)
    col 17->20   : Global timestamp (in 25ns counts)
    col 21->24   : BxCounter (since the last resync)
    
    Line is there only if one of the comparison shows up an error for a given BX. Therefore, for all the other BX we can use the expected word as a baseline for reference and CIC frames.
    '''

    if line.find('|')!=-1:  # Parsing the error line from SEU log file
        nline=line.replace('| ','')
        nline=nline.replace('|','')
        nline=nline.replace('[','')
        nline=nline.replace(']','')
        nline=nline.replace('\'','')
        nline=nline.replace(',','')
        lData=nline.split(' ')
        
        error_type=hex( int(lData[7],16) )[2:].rjust(4,'0')
        time_stamp_1=int(lData[19]+lData[18]+lData[17]+lData[16],16)
        time_stamp_2=int(lData[23]+lData[22]+lData[21]+lData[20],16)
        
        treal=time_stamp_1*25*1e-9
        
        #print(treal)
        tcor=(time_stamp_2-19)%(ncycles*3564) # Frame in the firmware is in advance wrt the CIC one
        # Delay here is 19 clock cycles but it can change by +/- 1
        # Basically time_stamp_1 never changes, whereas the other depends a lot on the RESYNC signals
        
        # Do we have a trigger error ?
        if error_type[0:2]=='01' or error_type[0:2]=='03' or error_type[0:2]=='07':

            #print(time_stamp_1,int((time_stamp_2)/(ncycles*3564)),error_type[0:2])

            # If so we compare the data recorded at the ouptut of the CIC with
            # the expectation. Reminder, when an error is recorded, the expected
            # work is also available
            line_r=[]  # Reference (Should be equivalent to TRGO data)
            line_d=[]  # Data from CIC

            line_rb=[]
            line_db=[]
            
            line_rbb=[]
            line_dbb=[]


            # Data arrangement
            #    012345
            #    ABCDEF
            #    GHIJKL
            #    MNOPQR
            #    ...
            #
            #


            # We first put the data into the same format than our reference frame
            # These input blocks contain 16 bits / per line / BX: AGM... for line 0 for example
            for jj in range(6): # Loop over output lines
                line_r.append(bin(int(lData[8+jj],16))[2:].rjust(16,'0'))
                line_d.append(bin(int(lData[jj],16))[2:].rjust(16,'0'))
                
            for jj in range(int(8*freq/320)): # Comparing 8/16 6-bits words, depending on the output freq
                wd=[] # Data line at 320/640M (6 bits)
                wr=[] # Ref line at 320/640M (6 bits)
                wdd=[]
                wrr=[]
                
                # The bits are rearranged per clock cycle
                if freq==320: # 320M case, the bits are doubled
                    for kk in range(6):
                        wr += line_r[kk][2*jj]
                        wd += line_d[kk][2*jj]
                        wrr += line_r[kk][2*jj+1]
                        wdd += line_d[kk][2*jj+1]
                else:
                    for kk in range(6):
                        wr += line_r[kk][jj]
                        wd += line_d[kk][jj]

                # line_[0] contains ABCDEF now
                # We have 8 lines at 320M and 16 lines at 640M

                line_rb.append("".join(wr))
                line_db.append("".join(wd))
                line_rbb.append("".join(wrr))
                line_dbb.append("".join(wdd))
                
                u=zip(line_rb[jj],line_db[jj]) # Compare the 2 lines
                count=0
                for i,j in u:
                    if i!=j:
                        bitflip_trg+=1
                        count+=1
                if verbose==1 and count>0:
                    outlog.write(str(line_rb[jj])+'/'+str(line_db[jj])+'/'+str(int((time_stamp_2)/(ncycles*3564)))+'/'+str((time_stamp_2)%(ncycles*3564))+'\n')
                #else:
                #    outlog.write(str(line_rb[jj])+'\n')
                check=0
                
                # Thoses lines are just a control at 320M, there is no data here at 640M
                u=zip(line_db[jj],line_dbb[jj])
                for i,j in u:
                    if i!=j:
                        check+=1
                        
                # And store the error for subsequent analysis
                # error contains 11 parameters:
                #
                # error[0] : input frames is repeated many times in the CIC, here we give the
                #            repetition ID.
                # error[1] : the BX ID within the input test vector
                # error[2] : the clock tick within the BXID (bet 0 and 8/16 depending on output freq)
                # error[3/4]: the 6 bits of the corresponding clk tick in ref/data respectively
                # error[5] : number of bit flips in the line
                # error[6] : error type
                # error[7] : not used
                # error[8/9/10] : debug params to detect output data transmission issues (phase aligment problems), only possible at 320 where we duplicate bits
                
                error=[int((time_stamp_2)/(ncycles*3564)),(time_stamp_2)%(ncycles*3564),jj,line_rb[jj],line_db[jj],count,error_type[0:2],(time_stamp_2)%(8),line_rbb[jj],line_dbb[jj],check,treal]
                errlist_trg.append(error)

#
# 2. Parsing is completed, one can now analyze the different errors and compare the
#    bad frames to the reference ones
#

#print(errlist_trg)

errbycycles_trg=[]

# We first arrange the data by cycle
# More than one error in one pass of the input stream (a cycle) is possibly suspect
#

for error in errlist_trg:
    cycle=error[0]
    found=0
    #print(cycle,error)
    for data in errbycycles_trg:
        if cycle==data[0]:
            data.append(error)
            found=1
            break
    if found==0: # New cycle
        errbycycles_trg.append([cycle,error])

nflip_trg=0
nhead_error=0
nwidth_flip=0
naddedstubs=0
notinref=0
notindat=0

simple_err=0

simpleflip=0
simpleflip_pad=0

outlog.write('Number of LHC cycles in the input stream: '+str(ncycles)+'\n')

for cycle in errbycycles_trg:

    if verbose==1:
        outlog.write(str(cycle)+'\n')
    

    if len(cycle)-1>1: # Let's look at everything

        lines=0
        outlog.write('\n')
        outlog.write('//////////////////////////////////'+'\n')
        outlog.write('Start detailed analysis for cycle '+str(cycle[0])+' / '+str(cycle[1][1]%3564)+'\n')
        outlog.write('T_error = '+str(cycle[1][11])+'s\n')
        if len(cycle)-1==8*freq/320: # The simplest case (usually a flip in the output)
            outlog.write('Only one BX in error in this cycle'+'\n')
            simple_err+=1
            # Many single bitflip errors in different BXs in one cycle are highly unlikely
            # On the contrary it's likely to happen in case of IR drop issue
        outlog.write('\n')
        
        #
        # Here we make a very detailed comparison of the corresponding CIC
        # output block. We retrieve the expected output and store if in ref
        # and for the CIC output, when there is no error we take the expected
        # output, and the CIC one when there is a mismatch
        # This way we can reconstruct the complete data block, and check whether
        # the error is coming from an input data flip, or anything else
        #
        
        ref='' # The corresponding reference TRG output word
        dat='' # The CIC TRG output
        dat2=''# The same output, for debug purpose at 320M (phase alignment)
        # A new CIC output TRG block starts every 8BX, so we need to catch it here
        start_fr=cycle[1][1]-cycle[1][1]%8 # The start of the output block where the error is
        lgth_fr=int(max(64*freq/320,8*freq/320*(cycle[len(cycle)-1][1]-start_fr+1))) # Clk counts at 320/640
        if verbose==1:
            outlog.write('Properties of the output block retrieved (BXID/NLines/NBX): '+str(start_fr)+'/'+str(lgth_fr)+'/'+str(len(cycle)-1)+'\n')
        
        idx=[]   # Where to look in the TRGO vector
        false=[] # Will be used to store the wrong words ID
        check=0  # To detect output data transmission issues (only possible at 320M)

        
        # First of all we retrieve the reference data stream (stored in the bin file)
        
        for i in range(lgth_fr):
            ref=ref+TRGO[int(8*freq/320*start_fr+i)]
            idx.append(int(8*freq/320*start_fr+i))
            false.append(0)
            
        # Ref word is written, we now flag the blocks which are false
        # to correctly assign them in the data word
            
        for i in range(len(cycle)-1):
            if cycle[i+1][10]>0: # Potential output phase alignment issue for this block
                check=1          # require a special analysis
            for j in range(len(idx)):
                if idx[j]==int(8*freq/320*cycle[i+1][1]+cycle[i+1][2]):
                    false[j]=1
                    break
                    
        
        # We are now ready to write the data word
        
        nwrg=0

        if verbose==1:
            outlog.write('EXPECT//MEASURE (if different)'+'\n')

        for i in range(lgth_fr):

            if false[i]==1:
                dat=dat+cycle[nwrg+1][4]
                dat2=dat2+cycle[nwrg+1][9]
                if verbose==1:
                    outlog.write(str(cycle[nwrg+1][4])+'//'+str(TRGO[int(8*freq/320*start_fr+i)])+'\n')
                nwrg+=1
            else:
                if verbose==1:
                    outlog.write(str(TRGO[int(8*freq/320*start_fr+i)])+'\n')
                dat=dat+TRGO[int(8*freq/320*start_fr+i)]
                dat2=dat2+TRGO[int(8*freq/320*start_fr+i)]

        # Here we write the corresponding TRG output blocks retrieved
        
        for i in range(1):
            outlog.write('REF :'+str(ref[int(384*freq/320*i):int(384*freq/320*(i+1))])+'\n')
            outlog.write('DATA:'+str(dat[int(384*freq/320*i):int(384*freq/320*(i+1))])+'\n')
            if check!=0:
                outlog.write('Output phase issue '+str(dat2[int(384*freq/320*i):int(384*freq/320*(i+1))])+'\n')
            outlog.write('\n')
        
        
        #
        # Start of the detailed block analysis
        #
        
        # We analyze the number of stubs, and compare ref to data stubs
        # Bit flips could indeed come from a simple stub reordering in the output
        # word
        
        lstub=18   # Stub size, in bits (CBC case)
        if MPA==1:
            lstub=21

        nstubs=int((384*freq/320-28)/lstub) # Max number of stubs in the block

        refcnt=[]  # Stubs retrieved in the ref block
        datcnt=[]  # Stubs retrieved in the data block

        if verbose==1:
            outlog.write(str(cycle)+'\n')
            outlog.write(str(ref)+'\n')
            outlog.write(str(dat)+'\n')
            
        for i in range(nstubs): # Reference stubs
            nbits=0
            for p in ref[28+lstub*i:28+lstub*(i+1)]:
                if int(p)==1:
                    nbits+=1
            if nbits==0:  # There is no stub 00..00 so if we got there it's the end
                continue
            refcnt.append(ref[28+lstub*i:28+lstub*(i+1)])
            if verbose==1:
                outlog.write('Refstub '+str(i)+' '+str(ref[28+lstub*i:28+lstub*(i+1)])+'\n')
        outlog.write('\n')


        for i in range(nstubs): # Data stubs
            nbits=0
            for p in dat[28+lstub*i:28+lstub*(i+1)]:
                if int(p)==1:
                    nbits+=1
            if nbits==0:  # Here we can possibly have a bit flip, to be verified.
                continue
            datcnt.append(dat[28+lstub*i:28+lstub*(i+1)])
            if verbose==1:
                outlog.write('Datstub '+str(i)+' '+str(dat[28+lstub*i:28+lstub*(i+1)])+'\n')
        outlog.write('\n')
        
        # Compute the number of stubs stored in the blocks
        refmult=int(ref[27])+int(ref[26])*2+int(ref[25])*4+int(ref[24])*8+int(ref[23])*16+int(ref[22])*32
        datmult=int(dat[27])+int(dat[26])*2+int(dat[25])*4+int(dat[24])*8+int(dat[23])*16+int(dat[22])*32
        
        overflow=0
        if refmult==nstubs: # Overflow case is important, because different stubs can be thrown out
            overflow=1      # of the frame (bitonic sort output is not reproducible)
                            # so there could be differences bet ref and dat in such case
                            # which are not due to SEU
                            

        # Check for errors in the header (>>>we don't like them...<<<)
        if ref[0:22]!=dat[0:22]:
            if ref[0]!=dat[0]:
                outlog.write("Header error (CHIP TYPE)"+'\n')
            if ref[1:10]!=dat[1:10]:
                outlog.write("Header error (Error flags)"+'\n')
                outlog.write("Check the number of flags raised"+'\n')
                outlog.write("If 1: a synchro bit has likely flip in the corresponding chip"+'\n')
                outlog.write("If 8: likely a clock issue in the input"+'\n')
                outlog.write("!!Other cases should be properly investigated!!"+'\n')
                outlog.write('\n')
            if ref[11:22]!=dat[11:22]:
                outlog.write("Header error (BXID)"+'\n')
            outlog.write("REF HEADER: "+str(ref[0:22])+'\n')
            outlog.write("CIC HEADER: "+str(dat[0:22])+'\n')
            lines=1
            nhead_error+=1
            for data in cycle:
                outlog.write(str(data)+'\n')

        # Error in stub multiplicity is annoying, but less worrying
        # as stub can be artificially created by a bit flip. So Ndat>Nref is possible
        # The opposite a bit more dodgy
        multdiff=0
        if ref[22:28]!=dat[22:28]:
            outlog.write("Stub mult difference"+'\n')
            outlog.write("REF MULT: "+str(ref[22:28])+' = '+str(refmult)+'\n')
            outlog.write("CIC MULT: "+str(dat[22:28])+' = '+str(datmult)+'\n')
            outlog.write("If CIC-REF=1 the problem likely comes from a stub created in the payload"+'\n')
            outlog.write("Other situations should be analyzed carefully"+'\n')
            outlog.write('\n')
            lines=1
            multdiff=datmult-refmult

        # Then compare the stubs of ref and dat words

        missing_stubs=[]  # Stub in ref but not in data
        new_stubs=[]      # Stub in data but not in ref

        # First loop over the reference stubs to check if the stub is in the output block
        for stub in refcnt:
            found=0
            for ostub in datcnt:
                if ostub==stub:
                    found=1
                    break
            if found==0:
                outlog.write(str(stub)+' is not present in data'+'\n')
                lines=1
                missing_stubs.append(stub)
                
        # Second loop, the other way round...
        for stub in datcnt:
            found=0
            for ostub in refcnt:
                if ostub==stub:
                    found=1
                    break
            if found==0: # We found a non expected stub
                nbits=0
                nactive=0
                for p in stub[0:lstub]: # Compute the number of bitflip in the total stub payload
                    if int(p)==1:
                        nactive+=1
                for p in stub[6:lstub]: # Compute the number of bitflip in the MPA/CBC level stub payload (wo BX and Chip IDs)
                    if int(p)==1:
                        nbits+=1
                # Check if this new stub can come from a simple bit flip in the
                # payload of one CBC/MPA
                #
                # Problem is described in the SEU note, an artificial stub which just one bit in its payload can arise. This is a very specific situation, logically much less likely in CIC2.1. In MPA, by definition it can occur only in the 2nd BX of the input block
                
                if int(stub[2:3])==1 and nactive>1 and nbits==1 and MPA==1 and ((multdiff==0 and overflow==1) or (multdiff!=0 and overflow==0)):
                    outlog.write('Stubs come from a bit flip in input MPA word'+'\n')
                    lines=1
                    nwidth_flip+=1

                elif nbits==1 and MPA==0 and nactive>1 and multdiff==1 and ((multdiff==0 and overflow==1) or (multdiff!=0 and overflow==0)):
                    outlog.write('Stubs come from a bit flip in input CBC word'+'\n')
                    lines=1
                    nwidth_flip+=1

                if overflow==1:
                    outlog.write('In overflow so one good stub might be popped out'+'\n')

                if nactive>1:
                    new_stubs.append(stub)
                else:
                    outlog.write('Simple bit flip in padding bits'+'\n')
                    lines=1
                    simpleflip+=1
                    simpleflip_pad+=1

        #
        # Time to make a final comparison
        #

        notindat+=len(missing_stubs)
        notinref+=len(new_stubs)
        missing_stubs_final=[]
        new_stubs_final=[]
        
        # Compare all the new and missing stubs, if there is only a 1 bit difference an SEU flip in the output register is a pretty likely explanation
        for stubn in new_stubs:
            found=0
            for stubm in missing_stubs:
                u=zip(stubn,stubm)
                count=0
                for i,j in u:
                    if i!=j:
                        count+=1
                if count==1:
                    outlog.write('Simple bitflip in output (expected/recorded): '+str(stubm)+' / '+str(stubn)+'\n')
                    lines=1
                    found=1
                    simpleflip+=1
                    break
            if found==0:
                outlog.write('Brand new stub : '+str(stubn)+'\n')
                lines=1
        if lines==0:
            outlog.write('Simple bit flip in padding bits'+'\n')
            simpleflip+=1
            simpleflip_pad+=1
                
    
    
    # Number of bitflip in the other case (standard)
    #for i in range(len(cycle)-1):
    #    nflip_trg+=cycle[i+1][5]



#
# 3. Analysis is completed, present the results
#

outlog.write('\n')
outlog.write('_______________________________________________________'+'\n')
outlog.write('--> LogFile     : '+str(file)+'\n')
outlog.write('--> Input file  : '+fname+'\n')
outlog.write('--> Run length  : '+str(runtime)+'\n')

outlog.write('\n')
outlog.write('%%%% TRG PATH RESULTS %%%%'+'\n')
outlog.write('\n')
outlog.write('Total number of cycles with errors : '+str(len(errbycycles_trg))+'\n')
outlog.write('among which single/multi word errs : '+str(simple_err)+'/'+str(len(errbycycles_trg)-simple_err)+'\n')
outlog.write('Total number of bitflips           : '+str(bitflip_trg)+'\n')
outlog.write('Among which '+str(nhead_error)+' located in headers'+'\n')
outlog.write('Total number of missing stubs                                 :'+str(notindat)+'\n')
outlog.write('Total number of unexpected stubs                              :'+str(notinref)+'\n')
outlog.write('Corrected number of bitflips in the output word (in/out stubs):'+str(simpleflip)+'('+str(simpleflip-simpleflip_pad)+'/'+str(simpleflip_pad)+')'+'\n')
outlog.write('Remaining number of missing stubs                             :'+str(notindat-(simpleflip-simpleflip_pad))+'\n')
outlog.write('Remaining of stubs unexpected stubs                           :'+str(notinref-(simpleflip-simpleflip_pad))+'\n')


outlog.write('Among which coming from flips in the input (artificial stubs) :'+str(nwidth_flip)+'\n')
outlog.write('Unexplained of stubs unexpected stubs                         :'+str(notinref-(simpleflip-simpleflip_pad)-nwidth_flip)+'\n')

fp.close()
