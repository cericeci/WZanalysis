#!/bin/python
import fileinput,sys,os,subprocess,re,time
from datetime import datetime, timedelta
from random import choice
from os import listdir
from os.path import isfile, join
import time

pwd = os.getcwd()

listOfFiles=[]
typesList=[]
MClist=["ZZ", "Zgamma", "WV", "WZZJets", "ZZZJets", "WWZJets", "WWWJets", "TTWJets", "TTZJets", "TTWWJets", "TTGJets", "WWGJets"]
#MClist=[]
OtherList=["pu_syst", "ele_scale_syst", "mu_scale_syst", "ele_SF", "mu_SF"]
#DDlist=["dataDriven"]
DDlist=["dataDriven_el", "dataDriven_mu"]
METlist=["met_elEn", "met_jetEn", "met_jetRes", "met_muEn", "met_tauEn", "met_unEn"]
#METlist=[]
#ktermList=["kterm_up", "kterm_down"]
ktermDict={'kterm_up':6,'kterm_down':4}
other=["JER", "JES"]
#submit = True
submit = False
sysList = open('sysrun2.def', 'r').read().split('\n')
#sysList = open('sysrunTEST.def', 'r').read().split('\n')
outputNom="sys_nominal.def"
listOfFiles.append(outputNom)
sysList = sysList[:-1]
newfileNom=open(outputNom,'w')
for line in sysList:
    if line[:1]=="#":
        continue
    if len(line)==0:
        continue
    types= line.split()
    typesList.append(types[0])
    ####
    newfileNom.write(line+'\n')
    outputUp="sys_"+types[0]+"_UP.def"
    outputDown="sys_"+types[0]+"_DOWN.def"
    if types[0] in MClist:
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
variables=["Njets", "LeadingJetPt", "Zpt"]
#variables=["LeadingJetPt"]
algorithm="Bayes"
#global command

for v in variables:
    for var in variations:
        for nameMET in METlist:
            response= "sysResults/response_"+nameMET+"_"+var+".root"
            dataFile="inputHistograms2/forVukoNew/binning06/data_driven.root"
            WZfile="inputHistograms2/forVukoNew/binning06/"+nameMET+"_"+var+"/WZ.root"
            backgroundListName ="backgroundList06"
            outputFile="sysResults/unfolding_"+nameMET+"_"+v+"_"+var+".root"
            systFileName="sys_nominal.def"
            command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
            print command        
            if submit:
                p=subprocess.call(command,shell=True)
                
    for name in ktermDict:
        response="sysResults/response_sys_nominal.root"
        dataFile="inputHistograms2/forVukoNew/binning06/data_driven.root"
        WZfile= "inputHistograms2/forVukoNew/binning06/WZ.root"
        backgroundListName ="backgroundList06"
        outputFile="sysResults/unfolding_"+name+"_"+v+".root"
        command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -k "+ str(ktermDict[name])
        print command        
        if submit:
            p=subprocess.call(command,shell=True)

    for name in typesList:
        outputFile = "sysResults/unfolding_"+ name +"_"+ v + ".root"
        if name in DDlist:
            response= "sysResults/response_sys_nominal.root"
            dataFile= "inputHistograms2/forVukoNew/binning06/dataDrivenSys/data_driven_"+name+".root"
            WZfile= "inputHistograms2/forVukoNew/binning06/WZ.root"
            backgroundListName ="backgroundList06"
            outputFile = "sysResults/unfolding_"+ name +"_"+ v + ".root"
            systFileName="nominal.def"
            command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
            print command
            if submit:
                p = subprocess.call(command,shell=True)	
        else:
            for var in variations:
                if name in MClist:
                    response= "sysResults/response_sys_nominal.root"
                    dataFile= "inputHistograms2/forVukoNew/binning06/data_driven.root"
                    WZfile= "inputHistograms2/forVukoNew/binning06/WZ.root"
                    backgroundListName ="backgroundList_syst"
                    outputFile = "sysResults/unfolding_"+ name +"_"+ v +"_"+var+ ".root"
                    systFileName="sys_"+name+"_"+var+".def"
                    command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
                    print command
                    if submit:
                        p=subprocess.call(command, shell=True)
                if name in OtherList:
                    response= "sysResults/response_sys_"+name+"_"+var+".root"
                    dataFile= "inputHistograms2/forVukoNew/binning06/data_driven.root"
                    WZfile= "inputHistograms2/forVukoNew/binning06/"+name+"_"+var+"/WZ.root"
                    backgroundListName ="backgroundList06"+name+"_"+var
                    outputFile = "sysResults/unfolding_"+ name +"_"+ v +"_"+var+ ".root"
                    systFileName="sys_"+name+"_"+var+".def"
                    systFileName="sys_nominal.def"
                    command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
                    print command
                    if submit:
                        p=subprocess.call(command, shell=True)
                
                if name in other:
                    response= "sysResults/response_sys_"+name+"_"+var+".root"
                    dataFile= "inputHistograms2/forVukoNew/binning06/data_driven.root"
                    WZfile= "inputHistograms2/forVukoNew/binning06/WZ.root"
                    backgroundListName ="backgroundList06"
                    outputFile = "sysResults/unfolding_"+ name +"_"+ v +"_"+var+ ".root"
                    systFileName="sys_"+name+"_"+var+".def"
                    systFileName="sys_nominal.def"
                    command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
                    if submit:
                        p = subprocess.call(command,shell=True)	

    response= "sysResults/response_sys_nominal.root"
    dataFile= "inputHistograms2/forVukoNew/binning06/data_driven.root"
    WZfile= "inputHistograms2/forVukoNew/binning06/WZ.root"
    backgroundListName ="backgroundList06"
    outputFile = "sysResults/unfolding_nominal_"+v+".root"
    systFileName="sys_nominal.def"
    command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
    print command
    if submit:
        p = subprocess.call(command,shell=True)	


    response= "sysResults/response_met_nominal.root"
    dataFile="inputHistograms2/forVukoNew/binning06/data_driven.root"
    WZfile="inputHistograms2/forVukoNew/binning06/met_nominal/WZ.root"
    backgroundListName ="backgroundList06"
    outputFile="sysResults/unfolding_met_nominal_"+v+".root"
    systFileName="nominal.def"
    command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
    print command
    if submit:
        p = subprocess.call(command,shell=True)	
