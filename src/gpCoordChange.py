#!/usr/bin/python2.7

## A program to convert a the coordinates in a 
## GenePred file by adding or subtracting the desired
## adjustment.

## Written by: Kord M. Kober, kord.kober@ucsf.edu

import os,sys,subprocess,time,datetime
import re,getopt,string,array,math,random
import pkg_resources

OPERATOR="subtract"
DEBUG=0

def loadGenePredH(gpFile):
    if DEBUG: print "[loadGenePredH]",gpFile
    gpH = {}
    gpFP = open(gpFile,"r")
    for gpS in gpFP:
        gpA=gpS.strip("\n").split("\t")
        gpA.pop(0)  # remove the bin column
        if DEBUG: print "[loadGenePredH]",gpA
        gpH['name']     =gpA[0]
        gpH['chrom']    =gpA[1]
        gpH['strand']   =gpA[2]
        gpH['txStart']  = int(gpA[3])
        gpH['txEnd']    = int(gpA[4])
        gpH['cdsStart'] = int(gpA[5])
        gpH['cdsEnd']   = int(gpA[6])
        gpH['exonCount']    = int(gpA[7])
        gpH['exonStarts']   = []
        for start in gpA[8].split(','):
            if not start=="":
                gpH['exonStarts'].append(int(start))
        gpH['exonEnds']     = []
        for end in gpA[9].split(','):
            if not end=="":
                gpH['exonEnds'].append(int(end))
    if len(gpH['exonStarts']) != len(gpH['exonEnds']):
        print "[loadGenePredH] count of exonStarts != exonEnds"
        sys.exit(129)
    if len(gpH['exonStarts']) != gpH['exonCount']:
        print "[loadGenePredH] count of exonStarts != exonCount"
        sys.exit(130)
    else:
        if DEBUG: print "[loadGenePredH] Found",gpH['exonCount'],"exons."
    global GENEPREDNAME
    GENEPREDNAME = gpH['name']
    if DEBUG: print "[loadGenePredH]",gpH
    return gpH

def convertGenePredH(gpH,op,amt):
    if DEBUG: print "[convertGenePredH] in:",op,amt,gpH

    if op == "subtract":
        gpH['txStart']  = gpH['txStart'] - amt
        gpH['txEnd']    = gpH['txEnd'] - amt
        gpH['cdsStart'] = gpH['cdsStart'] - amt
        gpH['cdsEnd']   = gpH['cdsEnd'] - amt

        exonStartsA=[]
        for e in gpH['exonStarts']:
            exonStartsA.append(str(e - amt))
        gpH['exonStarts'] = exonStartsA

        exonEndsA=[]
        for e in gpH['exonEnds']:
            exonEndsA.append(str(e - amt))
        gpH['exonEnds'] = exonEndsA

    if DEBUG: print "[convertGenePredH] out:",gpH
    return gpH

#### Give some usage information for this program #######################################
def usage(errorNum):
    a = """
usage: gpCoordChange.py [hD] -p <genePredFile>
"""
    print a
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    genePredFile=""

    try:
        opts, args = getopt.getopt(argv, "hDp:s:a:", \
            ["help","debug","genePredFile","subtract","add"])
    except getopt.GetoptError:
        usage(20)
    for opt, arg in opts:
        global OPERATOR
        global DEBUG
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-p", "--genePredFile"):
            genePredFile = arg
        if opt in ("-s", "--subtract"):
            OPERATOR = "subtract"
            amount = int(arg)
        if opt in ("-a", "--add"):
            OPERATOR = "add"
            amount = int(arg)
        elif opt in ("-D", "--debug"):
            DEBUG = 1
            print "[main] Debug set."

    if genePredFile == "":
        usage(34)

    if DEBUG: print "[main] Will",OPERATOR,amount,"bases from the gene model coordinates."

    gpH_old=loadGenePredH(genePredFile)
    if DEBUG: print "[main] Using gene pred",gpH_old['name'],"from",genePredFile

    gpH=convertGenePredH(gpH_old,OPERATOR,amount)

    ## write it back out
    lineA=[]
    lineA.append(str(0))
    lineA.append(str(gpH["name"]))
    lineA.append(str(gpH["chrom"]))
    lineA.append(str(gpH["strand"]))
    lineA.append(str(gpH["txStart"]))
    lineA.append(str(gpH["txEnd"]))
    lineA.append(str(gpH["cdsStart"]))
    lineA.append(str(gpH["cdsEnd"]))
    lineA.append(str(gpH["exonCount"]))
    exonStartsS=",".join(gpH["exonStarts"])
    lineA.append(str(exonStartsS)+",")
    exonEndsS=",".join(gpH["exonEnds"])
    lineA.append(str(exonEndsS)+",")

    print("\t".join(lineA))


#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])

