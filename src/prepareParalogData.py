#!/usr/bin/python2.7

## A program to provide data to paralogTest.py

## Written by: Samantha Danison, samanthadanison00@gmail.com

import os,sys,subprocess,time,datetime
import re,getopt,string,array,math,random
import pkg_resources

DEBUG=0

#### Give some usage information for this program #######################################
def usage(errorNum):
    a = """
usage: prepareParalogData.py [-hD] -a <alignInfo.txt> -s <siteInfo.txt> -g <gene name>
"""
    print a
    if DEBUG: print "errorNum =",errorNum
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    alignInfoFile=""
    siteInfoFile=""
    gemeName=""

    try:
        opts,args = getopt.getopt(argv, "hDa:s:g:", \
            ["help","debug","alignInfoFile","siteInfoFile","geneName"])
    except getopt.GetoptError:
        usage(20)
    for opt, arg in opts:
        global DEBUG
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-a", "--alignInfoFile"):
            alignInfoFile = arg
        if opt in ("-s", "--siteInfoFile"):
            siteInfoFile = arg
        if opt in("-g","--geneName"):
            geneName=str(arg)
        elif opt in ("-D", "--debug"):
            DEBUG = 1
            print "[main] Debug set."

    if alignInfoFile == "":
        usage(34)
    if siteInfoFile == "":
        usage(35)
    if geneName == "":
        usage(36)

    dataDict={}

    alignFP = open(alignInfoFile, "r")
    for lineA in alignFP:
        line=lineA.strip("\n").split("\t")
        spp = line[0]
        if spp !='taxa':
            dataDict[spp]={}
            dataDict[spp]['averageCoverage']=line[7]
            dataDict[spp]['biMorphCount']=line[10]
            dataDict[spp]['triMorphCount']=line[8]
            dataDict[spp]['binLowCount']=0
            dataDict[spp]['binMidCount']=0
            dataDict[spp]['binHighCount']=0
    
    if DEBUG: print "After align stats added, dataDict=",dataDict 

    siteFP=open(siteInfoFile,"r")
    for lineA in siteFP:
        line=lineA.strip("\n").split("\t")
        spp=line[0]
        if spp != 'taxa':
            aFreq=float(line[2])
            cFreq=float(line[3])
            tFreq=float(line[4])
            gFreq=float(line[5])
            if not dataDict.has_key(spp):
                dataDict[spp]={}
            if aFreq>0.125 and aFreq<0.374:
                dataDict[spp]['binLowCount']+=1
            if aFreq>0.375 and aFreq<0.625:
                dataDict[spp]['binMidCount']+=1
            if aFreq>0.626 and aFreq<0.875:
                dataDict[spp]['binHighCount']+=1
            if cFreq>0.125 and cFreq<0.374:
                dataDict[spp]['binLowCount']+=1
            if cFreq>0.375 and cFreq<0.625:
                dataDict[spp]['binMidCount']+=1
            if cFreq>0.626 and cFreq<0.875:
                dataDict[spp]['binHighCount']+=1
            if tFreq>0.125 and tFreq<0.374:
                dataDict[spp]['binLowCount']+=1
            if tFreq>0.375 and tFreq<0.625:
                dataDict[spp]['binMidCount']+=1
            if tFreq>0.626 and tFreq<0.875:
                dataDict[spp]['binHighCount']+=1
            if gFreq>0.125 and gFreq<0.374:
                dataDict[spp]['binLowCount']+=1
            if gFreq>0.375 and gFreq<0.625:
                dataDict[spp]['binMidCount']+=1
            if gFreq>0.626 and gFreq<0.875:
                dataDict[spp]['binHighCount']+=1

    if DEBUG: print "After site info added dataDict=",dataDict

    headerA=['gene','spp','averageCoverage','biMorphCount','triMorphCount','binLowCount','binMidCount','binHighCount']
    paralogDataFP=open("alignStatsSummary.csv","w")
    paralogDataFP.write(",".join(headerA))
    paralogDataFP.write("\n")
    for spp,astatsH in dataDict.iteritems():
        lineA=[]
        lineA.append(geneName)
        lineA.append(spp)
        lineA.append(str(astatsH['averageCoverage']))
        lineA.append(str(astatsH['biMorphCount']))
        lineA.append(str(astatsH['triMorphCount']))
        lineA.append(str(astatsH['binLowCount']))
        lineA.append(str(astatsH['binMidCount']))
        lineA.append(str(astatsH['binHighCount']))
        paralogDataFP.write(",".join(lineA))
        paralogDataFP.write("\n")
    print "Wrote paralog data to file alignStatsSummary.csv"


#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])
