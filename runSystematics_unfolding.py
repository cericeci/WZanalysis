#!/bin/python
import fileinput,sys,os,subprocess,re,time
from datetime import datetime, timedelta
from random import choice
from os import listdir
from os.path import isfile, join

pwd = os.getcwd()

listOfFiles=[]
typesList=[]
MClist=["ZZ", "Zgamma", "WV", "WZZJets", "ZZZJets", "WWZJets", "WWWJets", "TTWJets", "TTZJets", "TTWWJets", "TTGJets", "WWGJets"]
OtherList=["pu_syst", "ele_scale_syst", "mu_scale_syst", "ele_SF", "mu_SF"]
DDlist=["dataDriven"]
other=["JER", "JES"]
submit = False
sysList = open('sysrun2.def', 'r').read().split('\n')
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

variations=["UP", "DOWN"]
variables=["Njets", "LeadingJetPt", "Zpt"]
algorithm="Bayes"
#global command

for v in variables:
    for name in typesList:
        outputFile = "sysResults/unfolding_"+ name +"_"+ v + ".root"
        if name in DDlist:
            response= "sysResults/response_nominal.root"
            dataFile= "inputHistograms2/forVukoNew/binning06/dataDrivenSyst/data_driven.root"
            WZfile= "inputHistograms2/forVukoNew/binning06/WZ.root"
            backgroundListName ="backgroundList06"
            systFileName="nominal.def"
            command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
            jobFileName = "jobs/unfolding-"+name+"_"+v+".csh"
            batchFileName = "jobs/unfolding-"+name+"_"+v+".bat"
            jobFile = open(jobFileName,"w")
            jobFile.write("#!/bin/tcsh \n")
            jobFile.write("cd " + pwd + "\n")
            jobFile.write(command + "\n");
            jobFile.close()

            batchFile = open(batchFileName,"w")
            batchFile.write("executable  =  "+jobFileName + '\n')
            batchFile.write("universe    =  vanilla \n")
            batchFile.write("log         =  condorlogs/unfolding-"+name+"_"+v+".log \n")
            batchFile.write("initialdir  =  "+pwd + "\n")
            batchFile.write("output      =  "+pwd+"/condorlogs/unfolding-"+name+"_"+v+".out \n")
            batchFile.write("error       =  "+pwd+"/condorlogs/unfolding-"+name+"_"+v+".err \n")
            batchFile.write("getenv      =  True \n\n")
            batchFile.write("queue \n")
            batchFile.close()
            
            p = subprocess.call("chmod +x " +jobFileName,shell=True)
            if submit:
                p = subprocess.call("condor_submit " +batchFileName,shell=True)	
        else:
            for var in variations:
                if name in MClist:
                    response= "sysResults/response_nominal.root"
                    dataFile= "inputHistograms2/forVukoNew/binning06/data_driven.root"
                    WZfile= "inputHistograms2/forVukoNew/binning06/WZ.root"
                    backgroundListName ="backgroundList"
                    outputFile = "sysResults/unfolding_"+ name +"_"+ v +"_"+var+ ".root"
                    systFileName="sys_"+name+"_"+var+".def"
                    command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
                #            print command
            
                if name in OtherList:
                    response= "sysResults/response_"+name+"_"+var+".root"
                    dataFile= "inputHistograms2/forVukoNew/binning06/data_driven.root"
                    WZfile= "inputHistograms2/forVukoNew/binning06/"+name+"_"+var+"/WZ.root"
                    backgroundListName ="backgroundList06"+name+"_"+var
                    outputFile = "sysResults/unfolding_"+ name +"_"+ v +"_"+var+ ".root"
                    systFileName="sys_"+name+"_"+var+".def"
                    command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName
                #print -command
                
                if name in other:
                    response= "sysResults/response_"+name+"_"+var+".root"
                    dataFile= "inputHistograms2/forVukoNew/binning06/data_driven.root"
                    WZfile= "inputHistograms2/forVukoNew/binning06/WZ.root"
                    backgroundListName ="backgroundList06"
                    outputFile = "sysResults/unfolding_"+ name +"_"+ v +"_"+var+ ".root"
                    #systFileName="sys_"+name+"_"+var+".def"
                    systFileName="sys_nominal.def"
                    command = "./wzDataUnfold -r "+ response + " -d " + dataFile + " -B "+ backgroundListName + " -t "+ WZfile+ " -v "+ v + " -a "+algorithm+" -o "+ outputFile + " -S "+ systFileName

                
                jobFileName = "jobs/unfolding-"+name+"_"+v+"_"+var+".csh"
                batchFileName = "jobs/unfolding-"+name+"_"+v+"_"+var+".bat"
                jobFile = open(jobFileName,"w")
                jobFile.write("#!/bin/tcsh \n")
                jobFile.write("cd " + pwd + "\n")
                jobFile.write(command + "\n");
                jobFile.close()

                batchFile = open(batchFileName,"w")
                batchFile.write("executable  =  "+jobFileName + '\n')
                batchFile.write("universe    =  vanilla \n")
                batchFile.write("log         =  condorlogs/unfolding-"+name+"_"+v+"_"+var+".log \n")
                batchFile.write("initialdir  =  "+pwd + "\n")
                batchFile.write("output      =  "+pwd+"/condorlogs/unfolding-"+name+"_"+v+"_"+var+".out \n")
                batchFile.write("error       =  "+pwd+"/condorlogs/unfolding-"+name+"_"+v+"_"+var+".err \n")
                batchFile.write("getenv      =  True \n\n")
                batchFile.write("queue \n")
                batchFile.close()
            
                p = subprocess.call("chmod +x " +jobFileName,shell=True)
                if submit:
                    p = subprocess.call("condor_submit " +batchFileName,shell=True)	
        
