'''

Data parser for CIC2 SEU test output file.
Takes as input a txt file containing the info for one run
Provides in output different type of info on TRG paths errors

'''



import sys, getopt
import os
import math
import pickle
import six

verbose=0

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:",["file="])
except getopt.GetoptError:
    print('DataParse.py -f <FileName>')
    sys.exit(2)

file=''

for opt, arg in opts:
    if opt == '-h':
        print('DataParse.py -f <FileName>')
        sys.exit()
    elif opt in ("-f", "--file"):
        file=arg

#
# 1. Open the file and parse the main info
#

freq=320

fp = open(file, 'r')

lines = fp.readlines()


errlist=[]
errlist_trg=[]

nerr=0
nherr=0
bitflip=0
MPA=1
t2=-1
runtime=0

ncycles=2
tcorprev=0
for line in lines:

    if line.find('elapsed')!=-1:
        runtime=(line.split(':')[2])
        
    if line.find('dumped')!=-1:
        print(line)
    if line.find('Analog')!=-1:
        print(line)

    if line.find('elapsed')!=-1:
        runtime=(line.split(':')[2])
    

    '''
    # Here we parse the error line provided by the analysis software:
    There are 24 16 bits words containing the following infos:
    For a bypass run only the four first trigger output lines are filled
    
    col 1->6     : CIC trigger lines output
    col 7        : L1 CIC output not realigned with RAM (unused)
    col 8        : Error flags (8 first bits for trigger, 8 following bits for L1).
    col 9->14    : RAM trigger lines (to be compared col 1->6)
    col 15       : RAM L1 line
    col 16       : CIC L1 output realigned (to be compared with 15)
    col 17->20   : Global timestamp (in 25ns counts)
    col 21->24   : BxCounter (since the last resync)
    
    Line is there only if one of the comparison shows up an error for a given BX. Therefore, for all the other BX we can use the expected word as a baseline for reference and CIC frames.
    '''
    firm_pb=0
    if line.find('|')!=-1:  # Parsing the error line from SEU log file
        nline=line.replace('| ','')
        nline=nline.replace('|','')
        nline=nline.replace('[','')
        nline=nline.replace(']','')
        nline=nline.replace('\'','')
        nline=nline.replace(',','')
        lData=nline.split(' ')
        

        #print(lData)
        error_type=hex( int(lData[7],16) )[2:].rjust(4,'0')
        time_stamp_1=int(lData[19]+lData[18]+lData[17]+lData[16],16)
        time_stamp_2=int(lData[23]+lData[22]+lData[21]+lData[20],16)
        
        tcor=(time_stamp_2-19)%(ncycles*3564) # Frame in the firmware is in advance wrt the CIC one
        #print lData
        # Bypass errors are considered as trigger error
        if error_type[0:2]=='01' or error_type[0:2]=='03' or error_type[0:2]=='07':

            line_r=[]  # Reference
            line_rrb=[]
            line_rb=[]
            line_d=[]  # Data from CIC
            line_db=[]
            line_ddb=[]
            
            # Bypassed data is output over the 4 first output lines
            for jj in range(4): # Loop over 4 first output lines
                line_r.append(bin(int(lData[8+jj],16))[2:].rjust(16,'0'))
                line_d.append(bin(int(lData[jj],16))[2:].rjust(16,'0'))
            
            for jj in range(4): # Comparing 8/16 6-bits words, depending on the output freq
                wd=[] # Data line at 320/640M (6 bits)
                wdd=[] # Data line at 320/640M (6 bits)
                wr=[] # Ref line at 320/640M (6 bits)
                wrr=[]
                
                firm_pb=0
                for kk in range(int(8*freq/320)):
                    wr  += line_r[jj][2*kk]
                    wrr += line_r[jj][2*kk+1]
                    wd  += line_d[jj][2*kk]
                    wdd += line_d[jj][2*kk+1]
                
                line_rb.append("".join(wr))
                line_rrb.append("".join(wrr))
                line_db.append("".join(wd))
                line_ddb.append("".join(wdd))

                u1=zip(line_rb[jj],line_db[jj])
                u2=zip(line_rb[jj],line_ddb[jj])
                u3=zip(line_rrb[jj],line_db[jj])
                u4=zip(line_rrb[jj],line_ddb[jj])
                
                count1=0
                count2=0
                count3=0
                count4=0
                mincount=10
                minidx=-1
                for i,j in u1:
                    if i!=j:
                        count1+=1

                if count1<mincount:
                    mincount=count1
                    minidx=1

                for i,j in u2:
                    if i!=j:
                        count2+=1
                      
                if count2<mincount:
                    mincount=count2
                    minidx=2
                      
                for i,j in u3:
                    if i!=j:
                        count3+=1
                      
                if count3<mincount:
                    mincount=count3
                    minidx=3
                    
                for i,j in u4:
                    if i!=j:
                        count4+=1
                      
                if count4<mincount:
                    mincount=count4
                    minidx=4
                
                # If at least one of the stream is reproducing the right input, then
                # the chip is OK
                #
                # Any difference at this level is related to an output phase misalignment in the
                # firmware, so we don't care
                
                if (line_db[jj]=='10101010' or line_db[jj]=='01010101' ):
                    continue
                if (line_ddb[jj]=='01010101' or line_ddb[jj]=='10101010'):
                    continue

                # If we are here it means we have a real error.
                # The tap thing are just there to realign the stream
                
                #print(tcorprev,tcor,jj,line_rb[jj],line_db[jj],line_ddb[jj],count)

                ur=zip(wr[1:8],wd[1:8])

                if mincount%2==1:
                    tap1="".join(wd[1:8])
                    tap2="".join(wdd[0:7])
                else:
                    tap1="".join(wd[0:7])
                    tap2="".join(wdd[1:8])
                    ur=zip(wr[1:8],wd[0:7])

                count_err=0
                
                for i,j in ur:
                    if i!=j:
                        count_err+=1

                print('')
                print("Error observed on I/O lines",jj," in cycle",int((time_stamp_2)/(ncycles*3564)),"at BX",(time_stamp_2)%(ncycles*3564))

                if (tap1==tap2):
                    if count_err==1:
                        nerr+=1
                    else:
                        print("Expected/Observed :","".join(wr[1:8]),"/",tap1)
                        print(count_err,"bitflips!!")
                        nherr+=1
                else:
                    print("Strange error",tap1,"/",tap2)
                    nherr+=1


#
#
# 2. Parsing is completed, one can now analyze the different errors and compare the
#    bad frames to the reference ones
#

print("Run duration in s",runtime)
print("Found",nerr,"simple error(s) corresponding to just one bitflip")
print("Found",nherr,"more complex error(s)")


errbycycles_trg=[]

nflip_trg=0
nstub_add_trg=0

# We first arrange the data by cycle

for error in errlist_trg:
    cycle=error[0]
    #    print error
    found=0
    for data in errbycycles_trg:
        if cycle==data[0]:
            data.append(error)
            found=1
            break
    if found==0:
        errbycycles_trg.append([cycle,error])

nhead_error=0
nwidth_flip=0
naddedstubs=0
notinref=0
notindat=0

simple_err=0

simpleflip=0
simpleflip_pad=0

sys.exit()
