#!/usr/bin/env python
import os, re
import commands
import math, time
import sys

print 
print 'START'
print 
print 'Submitting jobs for RadialMoment analysis '
print

# number of jobs to be submitted
NumberOfJobs=474
interval=4

FileList = "pPb_streamer_forests8TeV_2855##onwards.txt" # list with all the file directories

queue = "1nw" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
########   customization end   #########

path = os.getcwd()
print
print 'do not worry about folder creation:'
os.system("rm -r tmpdata")
os.system("mkdir tmpdata")
os.chdir("tmpdata")
print
for x in range(1, int(NumberOfJobs)+1):
    ##### creates directory and file list for job #######
    os.system("mkdir "+str(x))
    os.chdir(str(x))
    os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' ../../"+FileList+" > list.txt ")
    ##### creates jobs #######
    with open('job.sh', 'w') as fout:
        #fout.write("#!/bin/bash\n")
        fout.write("echo\n")
        fout.write("echo 'START---------------'\n")
        fout.write("echo 'WORKDIR ' ${PWD}\n")
        fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
        fout.write("pwd=`pwd`\n")
        fout.write("cd "+str(path)+"\n")
        fout.write("pwd\n")
        #fout.write("source setup.sh\n")
        fout.write("cp readForest_Data.exe tmpdata/"+str(x)+"/.\n")
        fout.write("cp Spring*DATA*.txt tmpdata/"+str(x)+"/.\n")
        fout.write("cd tmpdata/"+str(x)+"\n")
        fout.write("pwd\n")
        fout.write("cd /afs/cern.ch/work/r/rkunnawa/Run2_Analysis/RadialMoment/CMSSW_8_0_10/src\n")
        fout.write("echo 'Getting the proxy now'\n")
        fout.write("eval `scramv1 runtime -sh`\n")
        fout.write("export XRD_NETWORKSTACK=IPv4\n")
        fout.write("export X509_USER_PROXY=~/x509_user_proxy/proxy\n")
        fout.write("cd "+str(path)+"\n")
        fout.write("cd tmpdata/"+str(x)+"\n")
        fout.write("cat list.txt\n")
        fout.write("./readForest_Data.exe list.txt 0 "+str(interval)+" ak 4 PF 0 radialMoment_PythiaEpos_8TeVpthat50_part"+str(x)+".root\n")
        fout.write("echo 'STOP---------------'\n")
        fout.write("echo\n")
        fout.write("echo\n")
        
    os.system("chmod 744 job.sh")        
    ###### sends bjobs ######
    os.system("bsub -M 3000000 -q "+queue+" -o logs -J dataJob_"+str(x)+" < job.sh")
    print "job nr " + str(x) + " submitted"
        
    os.chdir("../")

print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print

