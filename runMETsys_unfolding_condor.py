#!/bin/python
import fileinput,sys,os,subprocess,re,time
from datetime import datetime, timedelta
from random import choice
from os import listdir
from os.path import isfile, join

pwd = os.getcwd()

submit = True
sysList = open('systMET', 'r').read().split('\n')

for l in sysList:
    name=l.split(".")

    if (name[0].split("_")[1]=="nominal"):
        outputFile="sysResults/response_met_"+name[0].split("_")[1]+".root"
        name_all=name[0].split("_")[1]
    else:
        if (name[0].split("_")[2]=="up"):
            outputFile="sysResults/response_met_"+name[0].split("_")[1]+"_UP.root"
            name_all=name[0].split("_")[1]+"_UP"
        if (name[0].split("_")[2]=="down"):
            outputFile="sysResults/response_met_"+name[0].split("_")[1]+"_DOWN.root"
            name_all= name[0].split("_")[1]+"_DOWN"
    command = "./wzMCUnfoldingAnalysis -l fullList_pucorr  -H Binnings/binning06 -o " + outputFile+" -S "+l
    print command
    print name_all
    jobFileName="jobs/response_met_"+name_all+".csh"
    batchFileName="jobs/response_met_"+name_all+".bat"
    jobFile = open(jobFileName,"w")
    jobFile.write("#!/bin/tcsh \n")
    jobFile.write("cd " + pwd + "\n")
    jobFile.write(command + "\n");
    jobFile.close()

    batchFile = open(batchFileName,"w")
    batchFile.write("executable  =  "+jobFileName + '\n')
    batchFile.write("universe    =  vanilla \n")
    batchFile.write("log         =  condorlogs/response-"+name_all+".log \n")
    batchFile.write("initialdir  =  "+pwd + "\n")
    batchFile.write("output      =  "+pwd+"/condorlogs/response-"+name_all+".out \n")
    batchFile.write("error       =  "+pwd+"/condorlogs/response-"+name_all+".err \n")
    batchFile.write("getenv      =  True \n\n")
    batchFile.write("queue \n")
    batchFile.close()

    p = subprocess.call("chmod +x " +jobFileName,shell=True)
    if submit:
        p = subprocess.call("condor_submit " +batchFileName,shell=True)	
 

