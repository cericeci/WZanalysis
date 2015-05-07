#!/bin/python
import fileinput,sys,os,subprocess,re,time
from datetime import datetime, timedelta
from random import choice
from os import listdir
from os.path import isfile, join
import math

nominalValuesList=open('filesForYieldsSystematics/numMC_nominal.h', 'r').read().split('\n')

muScaleDownList  =open('filesForYieldsSystematics/numMC_MUScale_DOWN.h', 'r').read().split('\n')
muScaleUpList    =open('filesForYieldsSystematics/numMC_MUScale_UP.h', 'r').read().split('\n')
elScaleDownList  =open('filesForYieldsSystematics/numMC_ELScale_DOWN.h', 'r').read().split('\n')
elScaleUpList    =open('filesForYieldsSystematics/numMC_ELScale_UP.h', 'r').read().split('\n')
muScaleFDownList  =open('filesForYieldsSystematics/numMC_MuScaleFactDOWN.h', 'r').read().split('\n')
muScaleFUpList    =open('filesForYieldsSystematics/numMC_MuScaleFactUP.h', 'r').read().split('\n')
elScaleFDownList  =open('filesForYieldsSystematics/numMC_ElScaleFactDOWN.h', 'r').read().split('\n')
elScaleFUpList    =open('filesForYieldsSystematics/numMC_ElScaleFactUP.h', 'r').read().split('\n')
pileupUpList      =open('filesForYieldsSystematics/numMC_pu_UP.h', 'r').read().split('\n')
pileupDownList      =open('filesForYieldsSystematics/numMC_pu_UP.h', 'r').read().split('\n')

ListUpSystematics  =[muScaleUpList,   elScaleUpList,   elScaleFUpList,   muScaleFUpList,   pileupUpList]
ListDownSystematics=[muScaleDownList, elScaleDownList, elScaleFDownList, muScaleFDownList, pileupDownList]
names=["muScale", "elScale", "ElScale factor", "MuScale factor", "pileup"]
totalSys={"ZZ3e":0,"ZZ2e1m":0,"ZZ1e2m":0,"ZZ3m":0,"Zg3e":0,"Zg2e1m":0,"Zg1e2m":0,"Zg3m":0,"WV3e":0,"WV2e1m":0,"WV1e2m":0,"WV3m":0,"VVV3e":0,"VVV2e1m":0,"VVV1e2m":0,"VVV3m":0}
#ZZ3e, ZZ2e1m, ZZ1e2m, ZZ3m, Zgamma3e, Zgamma2e1m, Zgamma1e2m, Zgamma3m, WV3e, WV2e1m, WV1e2m, WV3m, VVV3e, VVV2e1m, VVV1e2m, VVV3m
#for j in xrange(len(ListUpSystematics)):

for j in xrange(len(names)):
    print '***********************'
    print names[j]
    for i in xrange(len(nominalValuesList)):
        nominal=nominalValuesList[i].split(' ')
#        print "********"
#        print "i "+str(i)
#        print "j "+str(j)
        scaleDown=ListDownSystematics[j][i].split(' ')
        scaleUp  =ListUpSystematics[j][i].split(' ')

        if ((len(nominal)>2) and (len(scaleUp)>2) and (len(scaleDown)>2)):
            if (scaleUp[1][1]=="N"):
                if (float(nominal[2])<>0):
                    scaleSys=max( abs(float(scaleDown[2])-float(nominal[2])), abs(float(scaleUp[2])-float(nominal[2])))/float(nominal[2])
                    print nominal[1]+ "down: "+scaleDown[2]+" up: "+scaleUp[2]+" nominal: "+nominal[2]+" final: "+str(scaleSys)
                else:
                    scaleSys=0
                #print nominal[1]+'  '+str(scaleSys)#+ ' i:'+str(i)+' j '+str(j)+ ' k:'+str(k)+ ' ostatak: '+str(int(i+1)%16)
                if  (i==2):   
                    totalSys["ZZ3e"]=totalSys["ZZ3e"]+scaleSys*scaleSys
                if  (i==3):   
                    totalSys["ZZ2e1m"]=totalSys["ZZ2e1m"]+scaleSys*scaleSys
                if  (i==4):   
                    totalSys["ZZ1e2m"]=totalSys["ZZ1e2m"]+scaleSys*scaleSys
                if  (i==5):   
                    totalSys["ZZ3m"]=totalSys["ZZ3m"]+scaleSys*scaleSys
                if (i==10):
                    totalSys["Zg3e"]=totalSys["Zg3e"]+scaleSys*scaleSys
                if (i==11):
                    totalSys["Zg2e1m"]=totalSys["Zg2e1m"]+scaleSys*scaleSys
                if (i==12):
                    totalSys["Zg1e2m"]=totalSys["Zg1e2m"]+scaleSys*scaleSys
                    #print nominal[1]+'  '+str(scaleSys)#+ ' i:'+str(i)+' j '+str(j)+ ' k:'+str(k)+ ' ostatak: '+str(int(i+1)%16)
                if (i==13):
                    totalSys["Zg3m"]=totalSys["Zg3m"]+scaleSys*scaleSys
                if (i==18):
                    totalSys["WV3e"]=totalSys["WV3e"]+scaleSys*scaleSys
                if (i==19):
                    totalSys["WV2e1m"]=totalSys["WV2e1m"]+scaleSys*scaleSys
                if (i==20):
                    totalSys["WV1e2m"]=totalSys["WV1e2m"]+scaleSys*scaleSys
                if (i==21):
                    totalSys["WV3m"]=totalSys["WV3m"]+scaleSys*scaleSys
                    #print nominal[1]+'  '+str(scaleSys)#+ ' i:'+str(i)+' j '+str(j)+ ' k:'+str(k)+ ' ostatak: '+str(int(i+1)%16)
                if (i==26):
                    totalSys["VVV3e"]=totalSys["VVV3e"]+scaleSys*scaleSys
                if (i==27):
                    totalSys["VVV2e1m"]=totalSys["VVV2e1m"]+scaleSys*scaleSys
                if (i==28):
                    totalSys["VVV1e2m"]=totalSys["VVV1e2m"]+scaleSys*scaleSys
                if (i==29):
                    totalSys["VVV3m"]=totalSys["VVV3m"]+scaleSys*scaleSys
                
                #k=k+1
for key in totalSys:
#    print key
    print key+' '+str(math.sqrt(totalSys[key]))


