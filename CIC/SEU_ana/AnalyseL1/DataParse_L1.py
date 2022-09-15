'''

DataParse_L1.py

Main data parser for CIC2 SEU data (L1 part).

Takes as input a txt file containing the info for one run
Provides in output different type of info on L1 paths errors

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
import six

verbose=0  # Debug mode, put to 1 to get massive outlog.writeouts

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:",["file="])
except getopt.GetoptError:
    outlog.write('DataParse_L1.py -f <FileName>')
    sys.exit(2)

file=''

for opt, arg in opts:
    if opt == '-h':
        outlog.write('DataParse_L1.py -f <FileName>')
        sys.exit()
    elif opt in ("-f", "--file"):
        file=arg

#
# 1. Open the file and parse the main info
#

freq=640  # The default output freq is 640Mz

fsplit=file.split('/')
fname=fsplit[len(fsplit)-1]

outputlog='L1_SEU_ana_'+fname
outlog=open(outputlog,'w')

fp = open(file, 'r')

lines = fp.readlines()

lOutL1Bitstream=[]
llData=[]
L1O=[]
L1Oo=[]
L1Ostrip=[]

L1_starts=[]
sL1StartSeq = "1"*27+"0"
errlist=[]
errlist_trg=[]

bitflip_hard=0
bitflip=0
MPA=1
t2=-1
runtime=0

# Frame in the firmware is in advance wrt the CIC one
# delay is 20 for Louvain 2, 19 for Louvain 1
delay=0
shift=0
last=0
ncycles=4
hsize = 63 # Output L1 word header (MPA case)

#We loop over all the lines, and check for data

for line in lines:

    if line.find('elapsed')!=-1:
        runtime=(line.split(':')[2])
    
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
            hsize = 54
        if line.find('_640_')==-1:
            freq=320
        

        inputbare=fname.split("/")
        filename=inputbare[len(inputbare)-1]
        # filename is the name of the .bin file, in contains the main run
        # characteristics
        inputbare=filename.split(".")
        params=inputbare[0].split("_")

        ncycles = int(params[3])
        outlog.write('Opening bin file '+fname+'\n')



        inBitStreamFile = open(fname,'rb')
        input = pickle.load(inBitStreamFile) # Open the bit file w/pickle
        lOutL1Bitstream=input[3]             # The L1 reference output
        #llData=input[4]                      # The output stored on firmware (not used for the moment)

        # L1O contains the full triggger output over the ncycles generated cycles
        # Just one line so much simpler than in the trigger case

        for i in range(len(lOutL1Bitstream)):
            L1O.append(lOutL1Bitstream[i])
        #for i in range(len(llData)):
        #    L1Oo.append(llData[i])
            
        #print(L1Oo)
        sL1line = "".join(L1O) # The complete L1 output stream measured at the output of the CIC in normal case
        #sL1line_f = "".join(L1Oo) # The complete L1 output stream measured at the output of the CIC in normal case

        if verbose==1:
            outlog.write('Complete L1 output vector from the CIC '+sL1line+'\n')
            #outlog.write('Complete L1 output vector from the CIC '+sL1line_f+'\n')
        
        # L1 path is asynchronous, so we need to
        # look where the L1 blocks are starting
        # it's easy thanks to the unique start sequences
        # of the L1 output frame (27 bits at 1)
        
        lookstart=sL1line
        while sL1StartSeq in lookstart:
            startIndex = lookstart.find(sL1StartSeq)
            lookstart="0"*(startIndex+28)+lookstart[startIndex+28:]
            L1_starts.append(int(startIndex))

        if verbose==1:
            outlog.write('List of the starting BXIDs for L1 words '+str(L1_starts)+'\n')

        # Version of the L1 output formated like the firmware output
        # Hexa 16 bits words

        if freq==320: # For 320M the output bits are duplicated by the firmware
            for i in range(int(len(L1O)/8)):
                L1word=[]
                for j in range(8):
                    L1word.append(L1O[8*i+j])
                    L1word.append(L1O[8*i+j]) # Duplication (keep it to detect output phase alignment issues)
                sL1word="".join(L1word)
                L1Ostrip.append(hex( int(sL1word,2) )[2:].rjust(4,'0'))
        else:
            for i in range(int(len(L1O)/16)):
                L1word=[]
                for j in range(16):
                    L1word.append(L1O[16*i+j])
                sL1word="".join(L1word)
                L1Ostrip.append(hex( int(sL1word,2) )[2:].rjust(4,'0'))
             
        '''
        if verbose==1:
            outlog.write('L1 input vector (in hexa 16nits words) '+str(L1Ostrip)+'\n')
        '''
        
    # End of the input file opening loop
        
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
        
        tcor=(time_stamp_2-delay)%(ncycles*3564) # Frame in the firmware is in advance wrt the CIC one
        # delay is 20 for Louvain 2, 19 for Louvain 1

        ''' ???
        if tcor-last!=1:
            shift=0
    
        if lData[24]!='0000':
            shift+=1
            tcor=(time_stamp_2-delay-shift)%(ncycles*3564)
        '''
        # Do we have an L1 error ?
        if error_type[2:4]=='01' or error_type[2:4]=='03' or error_type[2:4]=='07':
            
            #if L1Ostrip[tcor]==lData[15]: # Known firmware bug
            #    continue                  # Should be corrected now but appear on old data

            # If so we compare the data recorded at the ouptut of the CIC with
            # the expectation. Reminder, when an error is recorded, the expected
            # CIC output is also available in L1Ostrip table

            rfwd=bin(int(L1Ostrip[tcor],16))[2:].rjust(16,'0') # Reference
            fiwd=bin(int(lData[14],16))[2:].rjust(16,'0')      # Firmware
            dfwd=bin(int(lData[15],16))[2:].rjust(16,'0')      # Data from CIC

            firword=''
            refword=''
            datword=''
            datword2=''
            phase=0
            #Put back the data line in 8/16bit binary format
            if freq==320:
                refword=rfwd[0]+rfwd[2]+rfwd[4]+rfwd[6]+rfwd[8]+rfwd[10]+rfwd[12]+rfwd[14]
                firword=fiwd[0]+fiwd[2]+fiwd[4]+fiwd[6]+fiwd[8]+fiwd[10]+fiwd[12]+fiwd[14]
                datword=dfwd[0]+dfwd[2]+dfwd[4]+dfwd[6]+dfwd[8]+dfwd[10]+dfwd[12]+dfwd[14]
                datword2=dfwd[1]+dfwd[3]+dfwd[5]+dfwd[7]+dfwd[9]+dfwd[11]+dfwd[13]+dfwd[15] # Control
            else:
                refword=rfwd[0]+rfwd[1]+rfwd[2]+rfwd[3]+rfwd[4]+rfwd[5]+rfwd[6]+rfwd[7]+rfwd[8]+rfwd[9]+rfwd[10]+rfwd[11]+rfwd[12]+rfwd[13]+rfwd[14]+rfwd[15]
                datword=dfwd[0]+dfwd[1]+dfwd[2]+dfwd[3]+dfwd[4]+dfwd[5]+dfwd[6]+dfwd[7]+dfwd[8]+dfwd[9]+dfwd[10]+dfwd[11]+dfwd[12]+dfwd[13]+dfwd[14]+dfwd[15]
                datword2=datword

            if datword!=datword2: # At 320M only both data should be identical otherwise it means one has a output line phase alignment problem (not CIC related a priori).
                phase=1
                outlog.write('Output phase alignment issue '+str(tcor)+' / '+str(shift)+' / '+str(refword)+' / '+str(datword)+' / '+str(datword2)+' / '+str(firword)+'\n')
            
            if verbose==1:
                outlog.write('Output params '+str(tcor)+' / '+str(shift)+' / '+str(refword)+' / '+str(datword)+' / '+str(datword2)+' / '+str(firword)+'\n')


            # Here we localize the error wrt the L1 words expected at the output

            countoi=8*freq/320*(tcor)  # countoi is the location of the error (in 320/640 clock ticks) in the original frame
            counto=countoi
            idx=0
            
            # Verify if the error is within the L1 word  or not
            # and define the params which will enable to relocate the error within the
            # original frame
            
            for i in range(len(L1_starts)):
                if i==0:
                    continue
                if countoi<L1_starts[i]: # The error between L1 word i-1 and L1 word i
                    counto=countoi-L1_starts[i-1] # How far it is from the start of L1 word i-1?
                    idx=i-1
                    break
            if counto==countoi: # Error is after the last L1 word
                counto=countoi-L1_starts[len(L1_starts)-1]


            # Then we compare the words

            u=zip(refword,datword)
            count=0 #
            l1ref=''
            l1read=''
            enough=0
            for i,j in u:
                if i!=j: # Bitflip
                    # If error is in the header we increment a special counter
                    # as this can be a problem
                    if counto>=0 and counto<hsize:
                        if verbose==1 and enough==0:
                            outlog.write('!!!!!!!!HEADER FLIP!!!!!!!!'+'\n')
                            outlog.write('Bitflip recorded '+str(counto)+' bit(s) after the start of the L1 word... '+str(idx)+' / '+ str(L1_starts[idx+1])+' / '+str(L1_starts[idx])+' / '+str(tcor)+'\n')
                            outlog.write('!!!!!!!!\n')
                        enough=1
                        bitflip_hard+=1
                    count+=1
                    bitflip+=1
                counto+=1
        
        
            # And store the error for subsequent analysis
            # error contains 12 parameters:
            #
            # error[0] : input frames is repeated many times in the CIC, here we give the
            #            repetition ID.
            # error[1] : the BX ID within the input test vector
            # error[2/3/9]: the reference/firmware/CIC 16bits word in error (hex format)
            # error[4,5,12] : the reference/CIC/firmware word (binary format)
            # error[6] : number of bit fliped
            # error[7] : Distance bet. the last flipped bit and the start of the L1 word
            # error[8] : Index of the L1 ID to which the flip belong (in the input frame)
            # error[10] : flag for flip in header
            # error[11] : flag for output phase alignment pb (320M only)
                
            error=[int((time_stamp_2-delay)/(ncycles*3564)),tcor,L1Ostrip[tcor],lData[15],refword,datword,count,counto,idx,lData[14],enough,phase,firword]
            errlist.append(error)


#
# 2. Parsing is completed, one can now analyze the different errors and compare the
#    bad frames to the reference ones
#

errbycycles=[]

# We first arrange the data by cycle
# More than one error in one pass of the input stream (a cycle) is possibly suspect
#

for error in errlist:
    cycle=error[0]
    found=0
    for data in errbycycles:
        if cycle==data[0]:
            data.append(error)
            found=1
            break
    if found==0:
        errbycycles.append([cycle,error])

# Do the analysis for the cycles containing more than 3 errors, everything if verbose

nflip=0
nheadrpb=0
nflip_pad=0

outlog.write('Number of LHC cycles in the input stream:'+str(ncycles)+'\n')

for cycle in errbycycles:

    if verbose==1:
        outlog.write(str(cycle)+'\n')
        
    hdr=0       # Problem in header
    bigmess=0   # Problem spanning over many data blocks
    
    # First of all flag the  issues
    if (cycle[1][7]>=0) and (cycle[1][7]<=hsize):
        nheadrpb+=1
        hdr=1
    
    if len(cycle)-1>10:
        bigmess=1


    # Other cases require a detailed analysis whatever the frame is
    outlog.write('\n')
    outlog.write('//////////////////////////////////'+'\n')
    outlog.write('Start detailed analysis for cycle '+str(cycle[0])+'\n')
    outlog.write('\n')

    if hdr==1:
        outlog.write('>>>>>> Error lies in the header'+'\n')
        outlog.write('\n')
        
    if bigmess==1 and hdr==0:
        outlog.write('>>>>>> \'Soft Overflow\' issue (aka non-standard issue)'+'\n')
        outlog.write('\n')
                
    #
    # Here we make a very detailed comparison of the corresponding CIC
    # output block. We retrieve the expected output and store if in ref
    # and for the CIC output, when there is no error we take the expected
    # output, and the CIC one when there is a mismatch
    # This way we can reconstruct the complete data block, and check whether
    # the error is coming from an input data flip, or anything else
    #
                
    ref=''  # A. The corresponding reference L1 output word
    dat=''  # B. The L1 output word out of the CIC
    firm='' # C. The L1 output word expected from the firmware (could differ from A)
        
    start_fr=cycle[1][1]-cycle[1][7]/(8*freq/320)-1 # The start of the output block where the error is
    lgth_fr=cycle[len(cycle)-1][1]-start_fr+60
    idx=[]  # Where to look in the L1Ostrip vector
    false=[]# Will be used to store the wrong words ID
        
    if cycle[1][6]!=1 or hdr!=0:
        if len(cycle)-1>=3 and verbose==1:
            for i in range(len(cycle)-1):
                outlog.write(str(i)+' | '+str(cycle[i+1][1])+' | '+str(cycle[i+1][4])+' | '+str(cycle[i+1][5])+' | '+str(cycle[i+1][12])+'\n')
        
        outlog.write('List of error messages:\n')
        outlog.write('i |  BX  | REFERENCE|   CIC    | FIRMWARE'+'\n')
        for i in range(len(cycle)-1):
            if cycle[i+1][4] != cycle[i+1][5]:
                outlog.write(str(i)+' | '+str(cycle[i+1][1])+' | '+str(cycle[i+1][4])+' | '+str(cycle[i+1][5])+' | '+str(cycle[i+1][12])+'\n')
        outlog.write('\n')
        
    # First of all we retrieve the reference data stream (stored in the bin file)

    for i in range(int(lgth_fr)):
        rfwd=bin(int(L1Ostrip[int(start_fr+i)],16))[2:].rjust(16,'0')
        refword=rfwd[0]+rfwd[2]+rfwd[4]+rfwd[6]+rfwd[8]+rfwd[10]+rfwd[12]+rfwd[14]
            
        if freq==640:
            refword=rfwd[0]+rfwd[1]+rfwd[2]+rfwd[3]+rfwd[4]+rfwd[5]+rfwd[6]+rfwd[7]+rfwd[8]+rfwd[9]+rfwd[10]+rfwd[11]+rfwd[12]+rfwd[13]+rfwd[14]+rfwd[15]

        ref=ref+refword
        idx.append(start_fr+i)
        false.append(0)
        if verbose==1:
            outlog.write(str(i)+' / '+str(refword)+'\n')
    
    # Ref word is written, we now flag the blocks which are false
    # to correctly assign them in the data word
    
    for i in range(len(cycle)-1):
        for j in range(len(idx)):
            if idx[j]==cycle[i+1][1]:
                false[j]=1
                break
                            
    # We are now ready to write the data word
    nwrg=0
    for i in range(int(lgth_fr)):

        if false[i]==1:
            dat=dat+cycle[nwrg+1][5]
            rfwd=bin(int(cycle[nwrg+1][9],16))[2:].rjust(16,'0')
            refword=rfwd[0]+rfwd[2]+rfwd[4]+rfwd[6]+rfwd[8]+rfwd[10]+rfwd[12]+rfwd[14]
            if freq==640:
                refword=rfwd[0]+rfwd[1]+rfwd[2]+rfwd[3]+rfwd[4]+rfwd[5]+rfwd[6]+rfwd[7]+rfwd[8]+rfwd[9]+rfwd[10]+rfwd[11]+rfwd[12]+rfwd[13]+rfwd[14]+rfwd[15]
            firm=firm+refword
            nwrg+=1
        else:
            rfwd=bin(int(L1Ostrip[int(start_fr+i)],16))[2:].rjust(16,'0')
            refword=rfwd[0]+rfwd[2]+rfwd[4]+rfwd[6]+rfwd[8]+rfwd[10]+rfwd[12]+rfwd[14]
            if freq==640:
                refword=rfwd[0]+rfwd[1]+rfwd[2]+rfwd[3]+rfwd[4]+rfwd[5]+rfwd[6]+rfwd[7]+rfwd[8]+rfwd[9]+rfwd[10]+rfwd[11]+rfwd[12]+rfwd[13]+rfwd[14]+rfwd[15]
            dat=dat+refword
            firm=firm+refword


    if cycle[1][7]<=28 and hdr==1:
        outlog.write('-> Flip(s) in the start sequence (the 28 first bits)'+'\n')
        hdr=2


    if cycle[1][7]>=29: # Do this analysis only if error is not in the start sequence
        r_s = ref.find(sL1StartSeq)
        d_s = dat.find(sL1StartSeq)

        p_size = 17  # P cluster size
        s_size = 14
            

        ref_pay = []
        dat_pay = []
        
        # First of all we compute the cluster multiplicities
        # Of course P cluster mult means something only for MPA....
            
        ref_mult_P=0
        ref_mult_S=int(ref[52+r_s])+int(ref[51+r_s])*2+int(ref[50+r_s])*4+int(ref[49+r_s])*8+int(ref[48+r_s])*16+int(ref[47+r_s])*32+int(ref[46+r_s])*64
        dat_mult_P=0
        dat_mult_S=int(dat[52+r_s])+int(dat[51+r_s])*2+int(dat[50+r_s])*4+int(dat[49+r_s])*8+int(dat[48+r_s])*16+int(dat[47+r_s])*32+int(dat[46+r_s])*64

        if MPA==1:
            ref_mult_P=int(ref[62+r_s])+int(ref[61+r_s])*2+int(ref[60+r_s])*4+int(ref[59+r_s])*8+int(ref[58+r_s])*16+int(ref[57+r_s])*32+int(ref[56+r_s])*64
            dat_mult_P=int(dat[62+r_s])+int(dat[61+r_s])*2+int(dat[60+r_s])*4+int(dat[59+r_s])*8+int(dat[58+r_s])*16+int(dat[57+r_s])*32+int(dat[56+r_s])*64

        # Then we retrieve the frames
        ref_pay=ref[r_s+hsize:r_s+hsize+p_size*ref_mult_P+s_size*ref_mult_S]
        dat_pay=dat[r_s+hsize:r_s+hsize+p_size*dat_mult_P+s_size*dat_mult_S]

        #print(str(ref_mult_S)+' / '+str(ref_mult_P)+' / '+str(cycle[0]))


        if ref[r_s+28:r_s+37]!=dat[r_s+28:r_s+37]:
            outlog.write('Status bit error'+'\n')
            hdr=2
        if ref[r_s+37:r_s+46]!=dat[r_s+37:r_s+46]:
            outlog.write('L1 ID error'+'\n')
            hdr=2


        if abs(ref_mult_S-dat_mult_S)>0:
            outlog.write('Strip clus mismatch (exp/rec): '+str(ref_mult_S)+'/'+str(dat_mult_S)+'\n')
        if abs(ref_mult_P-dat_mult_P)>0:
            outlog.write('Pixel clus mismatch (exp/rec): '+str(ref_mult_P)+'/'+str(dat_mult_P)+'\n')
                
        if verbose==1 or hdr>1:
            outlog.write('Reference stream  : '+str(ref)+'\n')
            outlog.write('CIC stream        : '+str(dat)+'\n')
            outlog.write('Firmware stream   : '+str(firm)+'\n')
            outlog.write('Reference header  : '+str(ref[r_s:r_s+hsize])+'\n')
            outlog.write('CIC header        : '+str(dat[r_s:r_s+hsize])+'\n')

            
        if cycle[1][6]!=1 or hdr!=0:
            outlog.write('Reference stream payload : '+str(ref_pay)+'\n')
            outlog.write('CIC stream payload       : '+str(dat_pay)+'\n')

        # Large cluster mult discrepancies are sign of a seriour problem, detailed analysis is needed then

        if abs(ref_mult_S-dat_mult_S)>2:
            outlog.write('Large S cluster mult discrepancy :'+str(ref_mult_S)+' / '+str(dat_mult_S)+'\n')
            outlog.write('Reference stream: '+str(ref)+'\n')
            outlog.write('CIC stream      : '+str(dat)+'\n')
            outlog.write('Firmware stream : '+str(firm)+'\n')
            continue

        if abs(ref_mult_P-dat_mult_P)>2:
            outlog.write('Large P cluster mult discrepancy :'+str(ref_mult_P)+' / '+str(dat_mult_P)+'\n')
            outlog.write('Reference stream: '+str(ref)+'\n')
            outlog.write('CIC stream      : '+str(dat)+'\n')
            outlog.write('Firmware stream : '+str(firm)+'\n')
            continue


        #
        # Start of the detailed block analysis
        #
        
        # We analyze the number of stubs, and compare ref to data stubs
        # Bit flips could indeed come from a simple stub reordering in the output
        # word

        refcnt=[] # Clusters retrieved in the ref block
        datcnt=[] # Clusters retrieved in the data block

        for i in range(ref_mult_S): # Reference Sclusts
            refcnt.append(ref_pay[s_size*i:s_size*(i+1)])

        for i in range(ref_mult_P): # Reference Pclusts
            refcnt.append(ref_pay[s_size*ref_mult_S+p_size*i:s_size*ref_mult_S+p_size*(i+1)])

        for i in range(dat_mult_S): # Data Sclusts
            datcnt.append(dat_pay[s_size*i:s_size*(i+1)])
        
        for i in range(dat_mult_P): # Data Pclusts
            datcnt.append(dat_pay[s_size*dat_mult_S+p_size*i:s_size*dat_mult_S+p_size*(i+1)])


        # Then compare the clusters of ref and dat words
            
        missing_cluss=[]  # Clusters in ref but not in data
        new_cluss=[]      # Clusters in data but not in ref

        # First loop over the reference clus to check if it is in the output block
        for clus in refcnt:
            found=0
            for oclus in datcnt:
                if oclus==clus:
                    found=1
                    break
            if found==0:
                #if verbose==1:
                outlog.write(str(clus)+' is not present in data'+'\n')
                missing_cluss.append(clus)
                
        # Second loop, the other way round...
        for clus in datcnt:
            found=0
            for oclus in refcnt:
                if oclus==clus:
                    found=1
                    break
            if found==0: # We found a non expected clus
                if verbose==1:
                    outlog.write(str(clus)+' is not expected'+'\n')
                new_cluss.append(clus)
                    
        #
        # Time to make a final comparison
        #
        
        lines=0
        missing_cluss_final=[]
        new_cluss_final=[]
            
        # Compare all the new and missing clus, if there is only a 1 bit difference an SEU flip in the output register is a pretty likely explanation
        for clusn in new_cluss:
            found=0
            for clusm in missing_cluss:
                u=zip(clusn,clusm)
                count=0
                for i,j in u:
                    if i!=j:
                        count+=1
                if count==1:
                    outlog.write('Simple bitflip in output (expected/recorded): '+str(clusm)+' / '+str(clusn)+'\n')
                    lines=1
                    found=1
                    nflip+=1
                    break
            if found==0:
                outlog.write('Brand new clus : '+str(clusn)+'\n')
                lines=1
        if lines==0:
            outlog.write('Simple bit flip in padding bits'+'\n')
            nflip+=1
            nflip_pad+=1



#
# 3. Analysis is completed, present the results
#

outlog.write('\n')
outlog.write('_______________________________________________________'+'\n')
outlog.write('--> LogFile     : '+str(file)+'\n')
outlog.write('--> Input file  : '+fname+'\n')
outlog.write('--> Run length  : '+str(runtime)+'\n')

outlog.write('\n')
outlog.write('%%%% L1 PATH RESULTS %%%%'+'\n')
outlog.write('\n')


outlog.write('Total number of bitflips '+str(bitflip)+'\n')
outlog.write('Among which '+str(bitflip_hard)+' located in headers'+'\n')
outlog.write('Corrected number of bitflips '+str(nflip)+'\n')
outlog.write('Among which '+str(nheadrpb)+' located in headers'+'\n')

fp.close()
