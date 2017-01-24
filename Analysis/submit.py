#!/usr/bin/env python
import os, re
import commands
import math, time
import sys

print 
print 'START'
print 

# total number of files = 8194: 34 * 241

# location of output in eos
outdir="root://eosuser.cern.ch://eos/user/r/rkunnawa/Run2_pPbData/"
# number of jobs to be submitted
NumberOfJobs= 10
# number files to be processed in a single job, take care to split your file so that you run on all files. The last job might be with smaller number of files (the ones that remain).
interval = 10
# base of the output file name, they will be saved in res directory
OutputFileNames = "pPbData_MinBias_8TeV_histograms"
radius=4
algo=""
jettype="PF"
# script to be used with cmsRun
FileList = "pPb_MinBias8TeV_forests.txt" # list with all the file directories
queue = "2nd" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
########   customization end   #########

path = os.getcwd()
print
print 'do not worry about folder creation:'
os.system("rm -r tmpdata")
os.system("mkdir tmpdata")
os.system("mkdir resdata")
print

##### loop for creating and sending jobs #####
for x in range(1, int(NumberOfJobs)+1):
   ##### creates directory and file list for job #######
   os.system("mkdir tmpdata/"+str(x))
   os.chdir("tmpdata/"+str(x))
   os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' ../../"+FileList+" > list.txt ")
   
   ##### creates jobs #######
   with open('job.sh', 'w') as fout:
      #fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo 'START---------------'\n")
      fout.write("echo 'WORKDIR ' ${PWD}\n")
      fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
      fout.write("pwd=`pwd`\n")
      fout.write("cd /afs/cern.ch/work/r/rkunnawa/Run2_Analysis/QUARKvsGLUON/ForestMaking/CMSSW_8_0_23/src\n")
      fout.write("pwd\n")
      fout.write("cmsenv\n")
      fout.write("eval `scramv1 runtime -sh`\n")
      fout.write("echo 'Getting the proxy now'\n")
      fout.write("export XRD_NETWORKSTACK=IPv4\n")
      fout.write("export X509_USER_PROXY=~/x509_user_proxy/proxy\n")
      fout.write("cd "+str(path)+"/tmpdata/"+str(x)+"\n")
      fout.write("cp /afs/cern.ch/work/r/rkunnawa/Run2_Analysis/RadialMoment/Analysis/makeHistograms_Data.C .\n")
      fout.write("root -b -l <<EOF\n")
      fout.write(".x makeHistograms_Data.C+("+str(1+interval*(x-1))+","+str(interval*x)+","+str(radius)+",\"" + algo + "\", \"" + jettype + "\",\"/afs/cern.ch/work/r/rkunnawa/Run2_Analysis/RadialMoment/Analysis/resdata/" + OutputFileNames+"_"+str(x)+".root" + "\")\n")
      fout.write(".q\n")
      fout.write("EOF\n")
      fout.write("echo 'STOP---------------'\n")
      fout.write("echo\n")
      fout.write("echo\n")
   os.system("chmod 744 job.sh")
   
   ###### sends bjobs ######
   os.system("bsub -M 4000000 -q "+queue+" -o logs -J dataJob_"+str(x)+" < job.sh")
   print "job nr " + str(x) + " submitted"
   
   os.chdir("../..")
   
print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print
