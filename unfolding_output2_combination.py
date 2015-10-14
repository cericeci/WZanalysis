 #!/bin/python
import fileinput,sys,os,subprocess,re,time
from datetime import datetime, timedelta
from random import choice
from os import listdir
from os.path import isfile, join

def seeExponent(number):
    afterDot=False
    if (float(number)<1):
        if "e" in number:
            index=number.find("e")
            #print (number[int(index+3)]) 
            return -int(number[int(index+3)])
        for j in range (0, len(number)):
            if (number[j]=='.'):
                #afterDot=True
                continue
            if (int(number[j])==0):
                continue
            if (int(number[j])<>0):
                #if afterDot:
                return int(-(j-1))
            #else:
             #   return int(-j)
    else:
        for i in range (0, len(number)):
            if (number[i]=='.'):
                return i-1
    return 10000

#variable="Zpt"
#variable="LeadingJetPt"
variable="Njets"
names=["","","","","","","","","","", "","","",""]

namesZpt=["0-20 GeV", "20-40 GeV", "40-60 GeV", "60-80 GeV", "80-100 GeV", "100-120 GeV", "120-140 GeV", "140-200 GeV", "200-300 GeV"]
namesLeadingJetPt=["30-60 GeV","60-100 GeV", "100-150 GeV", "150-250 GeV"]
namesNjets=["0 jets", "1 jet", "2 jets", "3 jets", "4 jets"]
if (variable=="Zpt"):
    names=namesZpt
    linesMain= open("outCombination-Zpt.txt", 'r').read().split('\n')    
if (variable=="Njets"):
    names=namesNjets
    linesMain= open("outCombination-Njets.txt", 'r').read().split('\n')
if (variable=="LeadingJetPt"):
    names=namesLeadingJetPt
    linesMain= open("outCombination-LeadingJetPt.txt", 'r').read().split('\n')


#linesError1=open("outError1.txt", 'r').read().split('\n')
#linesError2=open("outError2.txt", 'r').read().split('\n')
#linesError3=open("outError3.txt", 'r').read().split('\n')

for i in range(0,len(linesMain)):
    base=10
    main_values=linesMain[i].split('\t')
 #   error1_values=linesError1[i].split(' ')
 #   error2_values=linesError2[i].split(' ')
 #   error3_values=linesError3[i].split(' ')
    stringForPrint=""
    if (len(main_values)>1):
        mainValue1=float(main_values[1])
        error1_1=float(main_values[2])
        error2_1=float(main_values[3])
        error3_1=float(main_values[4])
        mainExp1=seeExponent(main_values[1])
        error1Exp1=seeExponent(main_values[2])
        error2Exp1=seeExponent(main_values[3])
        error3Exp1=seeExponent(main_values[4])

#        mainValue2=float(main_values[2])
#        error1_2=float(error1_values[2])
#        error2_2=float(error2_values[2])
#        error3_2=float(error3_values[2])
#        mainExp2=seeExponent(main_values[2])
#        error1Exp2=seeExponent(error1_values[2])
#        error2Exp2=seeExponent(error2_values[2])
#        error3Exp2=seeExponent(error3_values[2])

#        mainValue3=float(main_values[3])
#        error1_3=float(error1_values[3])
#        error2_3=float(error2_values[3])
#        error3_3=float(error3_values[3])
#        mainExp3=seeExponent(main_values[3])
#        error1Exp3=seeExponent(error1_values[3])
#        error2Exp3=seeExponent(error2_values[3])
#        error3Exp3=seeExponent(error3_values[3])

 #       mainValue4=float(main_values[4])
 #       error1_4=float(error1_values[4])
 #       error2_4=float(error2_values[4])
 #       error3_4=float(error3_values[4])
 #       mainExp4=seeExponent(main_values[4])
 #       error1Exp4=seeExponent(error1_values[4])
 #       error2Exp4=seeExponent(error2_values[4])
 #       error3Exp4=seeExponent(error3_values[4])

#        mainValue5=float(main_values[5])
#        error1_5=float(error1_values[5])
#        error2_5=float(error2_values[5])
#        error3_5=float(error3_values[5])
#        mainExp5=seeExponent(main_values[5])
#        error1Exp5=seeExponent(error1_values[5])
#        error2Exp5=seeExponent(error2_values[5])
#        error3Exp5=seeExponent(error3_values[5])

        #print str(error2Exp1-
        rounding1=1
        rounding2=1
        rounding3=1
        rounding4=1
        rounding5=1
        
        if (mainValue1<1):
            rounding1=abs(error3Exp1)-abs(mainExp1)
        #if (mainValue2<1):
         #   rounding2=abs(error3Exp2)-abs(mainExp2)
        #if (mainValue3<1):
        #    rounding3=abs(error3Exp3)-abs(mainExp3)
        #if (mainValue4<1):
        #    rounding4=abs(error3Exp4)-abs(mainExp4)
        #if (mainValue5<1):
        #    rounding5=abs(error3Exp5)-abs(mainExp5)

        stringForPrint=names[i]+" & "
        #print names[i]+" & ",
        if ((abs(mainExp1)==1) or (abs(mainExp1)==0)) and (variable=="Njets"):
            main=str(format(mainValue1, "."+str(rounding1+1)+"f"))
            string1=" $\pm$ "+str(format(error1_1,"."+str(rounding1+1)+"f")) 
            string2=" $\pm$ "+str(format(error2_1,"."+str(rounding1+1)+"f"))
            string3=" $\pm$ "+str(format(error3_1,"."+str(rounding1+1)+"f")) +"  \\\\ "
            total=main+string1+string2+string3
            stringForPrint=stringForPrint+total
        else:
            main= "("+str(format(mainValue1*base**abs(mainExp1),"."+str(rounding1)+"f"))
            string1=" $\pm$ "+str(format(error1_1*base**abs(mainExp1),"."+ str(rounding1)+"f"))  
            string2=" $\pm$ "+str(format(error2_1*base**abs(mainExp1),"."+str(rounding1)+"f"))
            string3=" $\pm$ "+str(format(error3_1*base**abs(mainExp1),"."+str(rounding1)+"f"))+")$\cdot 10 ^{"+str(mainExp1)+"}$ \\\\ "
            total= main+string1+string2+string3
            stringForPrint=stringForPrint+total

#        if ((abs(mainExp2)==1) or (abs(mainExp2)==0)) and (variable=="Njets"):
#            main= str(format(mainValue2, "."+str(rounding2+1)+"f"))
#            string1=" $\pm $ "+str(format(error1_2,"."+str(rounding2+1)+"f")) 
#            string2=" $\pm $ "+str(format(error2_2,"."+str(rounding2+1)+"f"))
#            string3=" $\pm $ "+str(format(error3_2,"."+str(rounding2+1)+"f"))+" & "
#            total= main+string1+string2+string3
#            stringForPrint=stringForPrint+total
#        else:
#            main= "("+str(format(mainValue2*base**abs(mainExp2),"."+str(rounding2)+"f"))
#            string1=" $\pm$ "+str(format(error1_2*base**abs(mainExp2),"."+str(rounding2)+"f"))
#            string2=" $\pm$ "+str(format(error2_2*base**abs(mainExp2),"."+str(rounding2)+"f"))
#            string3=" $\pm$ "+str(format(error3_2*base**abs(mainExp2),"."+str(rounding2)+"f"))+")$\cdot 10 ^{"+str(mainExp2)+"}$ & "
#            total= main+string1+string2+string3
#            stringForPrint=stringForPrint+total

#        if ((abs(mainExp3)==1) or (abs(mainExp3)==0)) and (variable=="Njets"):
#            main= str(format(mainValue3, "."+str(rounding3+1)+"f"))
#            string1=" $\pm $ "+str(format(error1_3, "."+str(rounding3+1)+"f")) 
#            string2=" $\pm $ "+str(format(error2_3,"."+str(rounding3+1)+"f"))
#            string3=" $\pm $ "+str(format(error3_3,"."+str(rounding3+1)+"f"))+" & "
#            total= main+string1+string2+string3
#            stringForPrint=stringForPrint+total
#        else:
#            main= "("+str(format(mainValue3*base**abs(mainExp3), "."+str(rounding3)+"f"))
#            string1=" $\pm$ "+str(format(error1_3*base**abs(mainExp3),"."+str(rounding3)+"f")) 
#            string2=" $\pm$ "+str(format(error2_3*base**abs(mainExp3),"."+str(rounding3)+"f"))
#            string3=" $\pm$ "+str(format(error3_3*base**abs(mainExp3),"."+str(rounding3)+"f"))+")$\cdot 10 ^{"+str(mainExp3)+"}$ & "
#            total= main+string1+string2+string3
#            stringForPrint=stringForPrint+total
#        if ((abs(mainExp4)==1) or (abs(mainExp4)==0)) and (variable=="Njets"):
#            main= str(format(mainValue4, "."+str(rounding4+1)+"f"))+" $\pm$ "
#            string1=" $\pm$ "+str(format(error1_4, "."+str(rounding4+1)+"f")) 
#            string2=" $\pm$ "+str(format(error2_4,"."+str(rounding4+1)+"f"))
#            string3=" $\pm$ "+str(format(error3_4,"."+str(rounding4+1)+"f"))+" & "
#            total=main+string1+string2+string3
#            stringForPrint=stringForPrint+total
#        else:
#            main= "("+str(format(mainValue4*base**abs(mainExp4), "."+str(rounding4)+"f"))
#            string1=" $\pm$ "+str(format(error1_4*base**abs(mainExp4),"."+str(rounding4)+"f")) 
#            string2=" $\pm$ "+str(format(error2_4*base**abs(mainExp4),"."+str(rounding4)+"f"))
#            string3=" $\pm$ "+str(format(error3_4*base**abs(mainExp4),"."+str(rounding4)+"f"))+")$\cdot 10 ^{"+str(mainExp4)+"}$ \\\\ "
#            total= main+string1+string2+string3
#            stringForPrint=stringForPrint+total

#        if (((abs(mainExp5)==1) or (abs(mainExp5)==0)) and (variable=="Njets")):
#            main = str(format(mainValue5, "."+str(rounding5+1)+"f"))
#            string1=" $\pm$ "+str(format(error1_5, "."+str(rounding5+1)+"f")) 
#            string2=" $\pm$ "+str(format(error2_5,"."+str(rounding5+1)+"f"))
#            string3=" $\pm$ "+str(format(error3_5,"."+str(rounding5+1)+"f"))+" \\\\"
#            total= main+string1+string2+string3
#            stringForPrint=stringForPrint+total
#        else:
#            main= "("+str(format(mainValue5*base**abs(mainExp5), "."+str(rounding5)+"f"))
#            string1=" $\pm$ "+str(format(error1_5*base**abs(mainExp5),"."+str(rounding5)+"f")) 
#            string2=" $\pm$ "+str(format(error2_5*base**abs(mainExp5),"."+str(rounding5)+"f"))
#            string3=" $\pm$ "+str(format(error3_5*base**abs(mainExp5),"."+str(rounding5)+"f"))+")$\cdot 10 ^{"+str(mainExp5)+"}$ \\\\"
#            total= main+string1+string2+string3
#            stringForPrint=stringForPrint+total
        print stringForPrint



#        if (mainExp4==1) or (mainExp4==0):
#            print str(format(mainValue4,"." +str(error2Exp4-mainExp4+1)+"f"))+" $\pm$ "+str(format(error1_4, "."+str(error2Exp4-mainExp4+1)+"f")) + "$\pm$ "+str(format(error2_4,"."+str(error2Exp4-mainExp4+1)+"f"))+" & ",
#        else:
#            print "("+str(format(mainValue4*base**mainExp4,"."+str(error2Exp4-mainExp4)+"f"))+" $\pm$ "+str(format(error1_4*base**mainExp4,"."+str(error2Exp4-mainExp4)+"f")) + "$\pm$ "+str(format(error2_4*base**mainExp4,"."+str(error2Exp4-mainExp4)+"f"))+")$\cdot 10 ^{-"+str(mainExp4)+"}$ & ",
#        if (mainExp5==1) or (mainExp5==0):
#            print str(format(mainValue5, "."+str(error2Exp5-mainExp5+1)+"f"))+" $\pm$ "+str(format(error1_5, "."+str(error2Exp5-mainExp5+1)+"f")) + "$\pm$ "+str(format(error2_5,"."+str(error2Exp5-mainExp5+1)+"f"))+"\\\\"
#        else:
#            print "("+str(format(mainValue5*base**mainExp5,"."+str(error2Exp5-mainExp5)+"f"))+" $\pm$ "+str(format(error1_5*base**mainExp5, "."+str(error2Exp5-mainExp5)+"f")) + "$\pm$ "+str(format(error2_5*base**mainExp5,"."+str(error2Exp5-mainExp5)+"f"))+")$\cdot 10 ^{-"+str(mainExp5)+"}$ \\\\" 
#        print "*************"


