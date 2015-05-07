#!/bin/python
import fileinput,sys,os,subprocess,re,time
from datetime import datetime, timedelta
from random import choice
from os import listdir
from os.path import isfile, join

pwd = os.getcwd()
listOfFiles=[]
typesList=[]
submit = True
sysList = open('sys.def', 'r').read().split('\n')
outputNom="sys_nominal.def"
#sysList=sysList[:-1]
newFileNom=open(outputNom, 'w')
for line in sysList:
    if line[:1]=="#":
        continue
    if len(line)==0:
        continue
    types=line.split()
    typesList.append(types[0])
    newFileNom.write(line+'\n')
    outputUp="sys_"+types[0]+"_UP.def"
    outputDown="sys_"+types[0]+"_DOWN.def"
    listOfFiles.append(outputDown)
    listOfFiles.append(outputUp)
    newfileUp=open(outputUp,'w')
    newfileDown=open(outputDown,'w')
    for line2 in sysList:
        if line2!=line:
            newfileUp.write(line2+'\n')
            newfileDown.write(line2+'\n')
        else:
            newfileUp.write(types[0]+'\t'+'1 \n')
            newfileDown.write(types[0]+'\t'+'-1 \n')
    newfileUp.close()
    newfileDown.close()

variations=["UP", "DOWN"]
for name in typesList:
    for var in variations:
        outputFile="/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/forVukoNew/binning06/"+name+"_"+var+"/"
        systFileName="sys_"+name+"_"+var+".def"
        command1="./wzAnalysisMC -S "+systFileName+" -f "+outputFile+" -H Binnings/binning06"
        command2="./wzAnalysisMCAll -S "+systFileName+" -f "+outputFile+ " -H Binnings/binning06"
        print command1
        print command2
        if submit:
            p=subprocess.call(command1, shell=True)
        if submit:
            p=subprocess.call(command2, shell=True)
#nominal
outputFile="/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/forVukoNew/binning06/"
systFileName="sys_nominal.def"
command1="./wzAnalysisMC -S "+systFileName+" -f "+outputFile+ " -H Binnings/binning06"
command2="./wzAnalysisMCAll -S "+systFileName+" -f "+outputFile+ " -H Binning/binning06"
print command1
print command2
if submit:
    p=subprocess.call(command1, shell=True)
if submit:
    p=subprocess.call(command2, shell=True)
    
