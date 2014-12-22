#!/bin/python
import fileinput,sys,os,subprocess,re,time
from datetime import datetime, timedelta
from random import choice
from os import listdir
from os.path import isfile, join

pwd = os.getcwd()

listOfFiles=[]
submit = True
sysList = open('sysrun.def', 'r').read().split('\n')
#sysList = open('sysrunTEST.def', 'r').read().split('\n')
outputNom="sys_nominal"+".def"
listOfFiles.append(outputNom)
sysList = sysList[:-1]
newfileNom=open(outputNom,'w')
for line in sysList:
    if line[:1]=="#":
        continue
    if len(line)==0:
        continue
    types= line.split()
    newfileNom.write(line+'\n')
    outputUp="sys_"+types[0]+"_UP.def"
    outputDown="sys_"+types[0]+"_DOWN.def"
#    print output


    listOfFiles.append(outputDown)
    listOfFiles.append(outputUp)
    

    newfileUp=open(outputUp,'w')
    newfileDown=open(outputDown,'w')
    for line2 in sysList:
        if line2!=line:
#            newfileNom.write(line2+'\n')
            newfileUp.write(line2+'\n')
            newfileDown.write(line2+'\n')
        else:
#            newfileNom.write(line2+'\n')
            newfileUp.write(types[0]+'\t'+'1 \n')
            newfileDown.write(types[0]+'\t'+'-1 \n')
#    for type in types:
#        print type+"\n"
#    print "****"
for l in listOfFiles:
    name=l.split('.')
    outputFile = "sysResults/response_"+ name[0] + ".root"
    command = "./wzMCUnfoldingAnalysis -l fullList_pucorr  -H Binnings/binning06 -o " + outputFile+" -S "+l
    print command
    if submit:
        p = subprocess.call(command,shell=True)	
