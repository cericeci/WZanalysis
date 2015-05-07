#!/bin/python
import fileinput,sys,os,subprocess,re,time
from datetime import datetime, timedelta
from random import choice
from os import listdir
from os.path import isfile, join
import math
pwd = os.getcwd()

nominal     =open('sysrun_nominal_met.h', 'r').read().split('\n')
elEnUp      =open('sysrun_elEn_up.h', 'r').read().split('\n')
elEnDown    =open('sysrun_elEn_down.h', 'r').read().split('\n')   
jetEnUp     =open('sysrun_jetEn_up.h', 'r').read().split('\n')
jetEnDown   =open('sysrun_jetEn_down.h', 'r').read().split('\n')   
jetResUp    =open('sysrun_jetRes_up.h', 'r').read().split('\n')
jetResDown  =open('sysrun_jetRes_down.h', 'r').read().split('\n')   
muEnUp      =open('sysrun_muEn_up.h', 'r').read().split('\n')
muEnDown    =open('sysrun_muEn_down.h', 'r').read().split('\n')   
tauEnUp     =open('sysrun_tauEn_up.h', 'r').read().split('\n')
tauEnDown   =open('sysrun_tauEn_down.h', 'r').read().split('\n')   
unclEnUp    =open('sysrun_unEn_up.h', 'r').read().split('\n')
unclEnDown  =open('sysrun_unEn_down.h', 'r').read().split('\n')   

listAxe=["dRatio3e", "dRatio2e1mu", "dRatio1e2mu", "dRatio3mu"];
nominalAxe=[0,0,0,0]
elEnUpAxe=[0,0,0,0]
elEnDownAxe=[0,0,0,0]
jetEnUpAxe=[0,0,0,0]
jetEnDownAxe=[0,0,0,0]
jetResUpAxe=[0,0,0,0]
jetResDownAxe=[0,0,0,0]
muEnUpAxe=[0,0,0,0]
muEnDownAxe=[0,0,0,0]
tauEnUpAxe=[0,0,0,0]
tauEnDownAxe=[0,0,0,0]
unclEnUpAxe=[0,0,0,0]
unclEnDownAxe=[0,0,0,0]

for i in xrange(len(nominal)-1):
    word= nominal[i].split(' ')
    word_elEnUp=elEnUp[i].split(' ')
    word_elEnDown=elEnDown[i].split(' ')
    word_jetEnUp= jetEnUp[i].split(' ')
    word_jetEnDown= jetEnDown[i].split(' ')
    word_jetResUp= jetResUp[i].split(' ')
    word_jetResDown =jetResDown[i].split(' ')
    word_muEnUp= muEnUp[i].split(' ')
    word_muEnDown= muEnDown[i].split(' ')
    word_tauEnUp= tauEnUp[i].split(' ')
    word_tauEnDown= tauEnDown[i].split(' ')
    word_unclEnUp= unclEnUp[i].split(' ')
    word_unclEnDown= unclEnDown[i].split(' ')
    if (len(word)<2):
        continue
    if word[1] in listAxe:
        if (word[1]=="dRatio3e"):
            nominalAxe[0]=float(word[2])
            elEnUpAxe[0] =float(word_elEnUp[2])
            elEnDownAxe[0]=float(word_elEnDown[2])
            jetEnUpAxe[0]=float(word_jetEnUp[2])
            jetEnDownAxe[0]=float(word_jetEnDown[2])
            jetResUpAxe[0]=float(word_jetResUp[2])
            jetResDownAxe[0]=float(word_jetResDown[2])
            muEnUpAxe[0]=float(word_muEnUp[2])
            muEnDownAxe[0]=float(word_muEnDown[2])
            tauEnUpAxe[0]=float(word_tauEnUp[2])
            tauEnDownAxe[0]=float(word_tauEnDown[2])
            unclEnUpAxe[0]=float(word_unclEnUp[2])
            unclEnDownAxe[0]=float(word_unclEnDown[2])
        if (word[1]=="dRatio2e1mu"):
            nominalAxe[1]=float(word[2])
            elEnUpAxe[1] =float(word_elEnUp[2])
            elEnDownAxe[1]=float(word_elEnDown[2])
            jetEnUpAxe[1]=float(word_jetEnUp[2])
            jetEnDownAxe[1]=float(word_jetEnDown[2])
            jetResUpAxe[1]=float(word_jetResUp[2])
            jetResDownAxe[1]=float(word_jetResDown[2])
            muEnUpAxe[1]=float(word_muEnUp[2])
            muEnDownAxe[1]=float(word_muEnDown[2])
            tauEnUpAxe[1]=float(word_tauEnUp[2])
            tauEnDownAxe[1]=float(word_tauEnDown[2])
            unclEnUpAxe[1]=float(word_unclEnUp[2])
            unclEnDownAxe[1]=float(word_unclEnDown[2])
        if (word[1]=="dRatio1e2mu"):
            nominalAxe[2]=float(word[2])
            elEnUpAxe[2] =float(word_elEnUp[2])
            elEnDownAxe[2]=float(word_elEnDown[2])
            jetEnUpAxe[2]=float(word_jetEnUp[2])
            jetEnDownAxe[2]=float(word_jetEnDown[2])
            jetResUpAxe[2]=float(word_jetResUp[2])
            jetResDownAxe[2]=float(word_jetResDown[2])
            muEnUpAxe[2]=float(word_muEnUp[2])
            muEnDownAxe[2]=float(word_muEnDown[2])
            tauEnUpAxe[2]=float(word_tauEnUp[2])
            tauEnDownAxe[2]=float(word_tauEnDown[2])
            unclEnUpAxe[2]=float(word_unclEnUp[2])
            unclEnDownAxe[2]=float(word_unclEnDown[2])
        if (word[1]=="dRatio3mu"):
            nominalAxe[3]=float(word[2])
            elEnUpAxe[3] =float(word_elEnUp[2])
            elEnDownAxe[3]=float(word_elEnDown[2])
            jetEnUpAxe[3]=float(word_jetEnUp[2])
            jetEnDownAxe[3]=float(word_jetEnDown[2])
            jetResUpAxe[3]=float(word_jetResUp[2])
            jetResDownAxe[3]=float(word_jetResDown[2])
            muEnUpAxe[3]=float(word_muEnUp[2])
            muEnDownAxe[3]=float(word_muEnDown[2])
            tauEnUpAxe[3]=float(word_tauEnUp[2])
            tauEnDownAxe[3]=float(word_tauEnDown[2])
            unclEnUpAxe[3]=float(word_unclEnUp[2])
            unclEnDownAxe[3]=float(word_unclEnDown[2])
print "*****"
print nominalAxe[0]
print nominalAxe[1]
print nominalAxe[2]
print nominalAxe[3]
print "********"
#comput syst
elEnSyst=[0,0,0,0]
jetEnSyst=[0,0,0,0]
jetResSyst=[0,0,0,0]
muEnSyst=[0,0,0,0]
tauEnSyst=[0,0,0,0]
unclEnSyst=[0,0,0,0]
all=[0,0,0,0]
for j in range(0,4):
    elEnSyst[j]=100* max(abs(elEnUpAxe[j]-nominalAxe[j]), abs(elEnDownAxe[j]-nominalAxe[j]))/nominalAxe[j]
    jetEnSyst[j] =100*max(abs(jetEnUpAxe[j]-nominalAxe[j]), abs(jetEnDownAxe[j]-nominalAxe[j]))/nominalAxe[j]
    jetResSyst[j]= 100*max(abs(jetResUpAxe[j]-nominalAxe[j]), abs(jetResDownAxe[j]-nominalAxe[j]))/nominalAxe[j]
    muEnSyst[j]=100*max(abs(muEnUpAxe[j]-nominalAxe[j]), abs(muEnDownAxe[j]- nominalAxe[j]))/nominalAxe[j]
    tauEnSyst[j]=100* max(abs(tauEnUpAxe[j]-nominalAxe[j]), abs(tauEnDownAxe[j]-nominalAxe[j]))/nominalAxe[j]
    unclEnSyst[j]=100* max(abs(unclEnUpAxe[j]-nominalAxe[j]), abs(unclEnDownAxe[j]-nominalAxe[j]))/nominalAxe[j]

    print "3e"
    print "Electron energy: "+str(format(elEnSyst[j],".2f"))
    print "Jet energy: "+str(jetEnSyst[j])
    print "Jet resolution: "+str(jetResSyst[j])
    print "Muon energy: "+str(muEnSyst[j])
    print "Tau energy: "+str(tauEnSyst[j])
    print "Unclustered energy: "+str(unclEnSyst[j])
    print "All: "+str(math.sqrt(elEnSyst[j]**2+jetEnSyst[j]**2+jetResSyst[j]**2+muEnSyst[j]**2+tauEnSyst[j]**2+unclEnSyst[j]**2)*100)
    all[j]=math.sqrt(elEnSyst[j]**2+jetEnSyst[j]**2+jetResSyst[j]**2+muEnSyst[j]**2+tauEnSyst[j]**2+unclEnSyst[j]**2)
latexOutput=True;

if (latexOutput):
    print "& 3e & 2e1mu & 1e2mu & 3mu \\\\"
    print "\hline"
    print "Electron energy & "+str(format(elEnSyst[0],".2f"))+" & "+str(format(elEnSyst[1],".2f"))+ " & "+str(format(elEnSyst[2],".2f"))+" & "+str(format(elEnSyst[3],".2f"))+"\\\\"
    print "Jet energy & "+str(format(jetEnSyst[0],".2f"))+" & "+str(format(jetEnSyst[1],".2f"))+ " & "+str(format(jetEnSyst[2],".2f"))+" & "+str(format(jetEnSyst[3],".2f"))+"\\\\"
    print "Jet resolution & "+str(format(jetResSyst[0],".2f"))+" & "+str(format(jetResSyst[1],".2f"))+ " & "+str(format(jetResSyst[2],".2f"))+" & "+str(format(jetResSyst[3],".2f"))+"\\\\"
    print "Muon energy & "+str(format(muEnSyst[0],".2f"))+" & "+str(format(muEnSyst[1],".2f"))+ " & "+str(format(muEnSyst[2],".2f"))+" & "+str(format(muEnSyst[3],".2f"))+"\\\\"
    print "Tau energy & "+str(format(tauEnSyst[0],".2f"))+" & "+str(format(tauEnSyst[1],".2f"))+ " & "+str(format(tauEnSyst[2],".2f"))+" & "+str(format(tauEnSyst[3],".2f"))+"\\\\"
    print "Unclustered energy & "+str(format(unclEnSyst[0],".2f"))+" & "+str(format(unclEnSyst[1],".2f"))+ " & "+str(format(unclEnSyst[2],".2f"))+" & "+str(format(unclEnSyst[3],".2f"))+"\\\\"
    print "\hline"
    print "All & "+str(format(all[0],".2f"))+" & "+str(format(all[1],".2f"))+ " & "+str(format(all[2],".2f"))+" & "+str(format(all[3],".2f"))+"\\\\"
    
