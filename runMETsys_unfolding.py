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
        command= "./wzAnalysisMC -S "+l+" -o "+name[0]
    else:
        if (name[0].split("_")[2]=="up"):
            outputFile="sysResults/response_met_"+name[0].split("_")[1]+"_UP.root"
        if (name[0].split("_")[2]=="down"):
            outputFile="sysResults/response_met_"+name[0].split("_")[1]+"_DOWN.root"
    command = "./wzMCUnfoldingAnalysis -l fullList_pucorr  -H Binnings/binning06 -o " + outputFile+" -S "+l
    print command
    if submit:
        p = subprocess.call(command,shell=True)	
