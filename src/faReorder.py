#!/usr/bin/python2.7

## A program to reorder a FASTA file given a set of labels

## Written by: Kord M. Kober, kord.kober@ucsf.edu

import os,sys,subprocess,time,datetime
import re,getopt,string,array,math,random
import pkg_resources
from Bio import SeqIO

DEBUG=0

#### Give some usage information for this program #######################################
def usage(errorNum):
    a = """
usage: faReorder.py [hD] -l <file with ordered labels> -i <input FASTA>
"""
    print a
    sys.exit(errorNum)

#### main #######################################
def main(argv):

    labelFileName=""
    inputFastaFileName=""

    try:
        opts, args = getopt.getopt(argv, "hDl:i:", \
            ["help","debug","labelFileName","inputFastaFileName"])
    except getopt.GetoptError:
        usage(20)
    for opt, arg in opts:
        global DEBUG
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-l", "--labelFileName"):
            labelFileName = arg
        if opt in ("-i", "--inputFastaFileName"):
            inputFastaFileName = arg
        elif opt in ("-D", "--debug"):
            DEBUG = 1
            print "[main] Debug set."

    if labelFileName == "":
        usage(34)
    if inputFastaFileName == "":
        usage(35)

    labelA=[]
    labelFP = open(labelFileName,"r")
    for labelS in labelFP:
        label=labelS.strip("\n")
        labelA.append(label)
    if DEBUG: print "[main] Loaded labels:",labelA

    faH={}
    for seq_record in SeqIO.parse(inputFastaFileName, "fasta"):
        if DEBUG: print "[main] Reading from FA:",seq_record.id
        faH[seq_record.id]=seq_record
    if DEBUG: print faH


    for l in labelA:    
        if faH.has_key(l):
            ## BioPython defaults to 60 character width
            print faH[l].format("fasta"),


#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])

