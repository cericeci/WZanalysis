#!/bin/python
import fileinput,sys,os,subprocess,re,time
from datetime import datetime, timedelta
from random import choice
from os import listdir
from os.path import isfile, join

pwd = os.getcwd()

#submit = False
submit = True
sysList = open('systMET', 'r').read().split('\n')

for l in sysList:
    name=l.split(".")

    if (name[0].split("_")[1]=="nominal"):
        outputPath="/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/forVukoNew/binning06/met_"+name[0].split("_")[1]+"/"
        #command= "./wzAnalysisMC -S "+l+" -o "+name[0]
    else:
        if (name[0].split("_")[2]=="up"):
            outputPath="/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/forVukoNew/binning06/met_"+name[0].split("_")[1]+"_UP/"
        if (name[0].split("_")[2]=="down"):
            outputPath="/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/forVukoNew/binning06/met_"+name[0].split("_")[1]+"_DOWN/"
    command= "./wzAnalysisMC -S "+l+" -o "+name[0]+".h -H Binning/binning06 -f "+outputPath
    print command
    if submit:
        p = subprocess.call(command,shell=True)	
