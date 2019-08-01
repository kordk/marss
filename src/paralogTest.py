#!/usr/bin/python2.7

import os,sys,subprocess,time,datetime
import re,getopt,string,array,math,random,numpy

#COVERAGE_FOLD_CUTOFF=4.0   ## depricated
BIN_DIFFERENCE_MID_MULTIPLIER=1
COVERAGE_STD_DEV_CUTOFF=2
COVERAGE_SPAL_DEPTH_CUTOFF=25.0

DEBUG=0

def coverageScore(testH,geneH):
    
    sppCovH={}
    #print geneH.keys()
    
    ## A. calculate the avgCovDepth and SD values using all genes for each spp.
    ##
    ## 1. collect all coverage values for each gene for each spp.
    for g,gStuffH in geneH.iteritems():
        #print "[coverageScore]",g,gStuffH.keys()
        for s,sStuffH in gStuffH.iteritems():
            #print "[coverageScore]",g,s,sStuffH['avgCovDepth']
            
            if not sppCovH.has_key(s):
                #print  "[coverageScore]",s
                sppCovH[s]={}
                sppCovH[s]['all']=[]
            sppCovH[s]['all'].append(sStuffH['avgCovDepth'])
            
    ## 2. determine mean and SD for each spp.
    for s,fooH in sppCovH.iteritems():
        sppCovH[s]['covAvg']=numpy.average(sppCovH[s]['all'])
        sppCovH[s]['covStd']=numpy.std(sppCovH[s]['all'])
        print  "[coverageScore]",s,"Avg(avgCovDepth)=",sppCovH[s]['covAvg'],"Std(avgCovDepth)=",sppCovH[s]['covStd']
    
    
    ## B. score each gene for coverage relative to spp. mean & SD
    numCoverageHighH={}
    for g,my_gStuffH in geneH.iteritems():
            
        for s,my_sStuffH in my_gStuffH.iteritems():
            
            if not testH[g].has_key(s):
                testH[g][s]={}
                testH[g][s]['coverageGt3SD']=0
                
            myDiffFromMean = abs(sppCovH[s]['covAvg'] - my_sStuffH['avgCovDepth'])
            if myDiffFromMean > sppCovH[s]['covStd']*COVERAGE_STD_DEV_CUTOFF:
                testH[g][s]['coverageGt3SD']+=1
                numCoverageHighH[g]=1
                
                if DEBUG: print "[coverageScore] High coverage observed:",g,s,"my avgCov=",my_sStuffH['avgCovDepth'],"spp. Mean=",sppCovH[s]['covAvg'],"spp. SD=",sppCovH[s]['covStd'],myDiffFromMean
            
            
        #sys.exit(44)
    
    print "[coverageScore] Found",len(numCoverageHighH),"genes where at least one species has an Avg(site avgCovDepth) >",COVERAGE_STD_DEV_CUTOFF,"SD"
    return testH,sppCovH


def trimorphicCounts(testH,geneH):
    
    numTrimorphicCountsH={}
    for g,gStuffH in geneH.iteritems():
            
        for s,sStuffH in gStuffH.iteritems():
            if not testH[g].has_key(s):
                testH[g][s]={}
            testH[g][s]['triMorphCount']=sStuffH['triMorphCount']
            
            #print "[trimorphicCounts]",g,s,testH[g][s]['triMorphCount']
            if testH[g][s]['triMorphCount'] > 0:
                numTrimorphicCountsH[g]=1
    
    print "[trimorphicCounts] Found", len(numTrimorphicCountsH),"genes with at least one species having at least one triMorphic site."
    return testH


def binCounts(testH,geneH):
    
    numSpalIgnoreH={}
    numInversedBinsAllH={}
    for g,gStuffH in geneH.iteritems():
        
        for s,sStuffH in gStuffH.iteritems():
            if not testH[g].has_key(s):
                testH[g][s]={}
                
            testH[g][s]['inversedBins']=0
                
            myLow=int(sStuffH['binLowCount'])
            myMid=int(sStuffH['binMidCount'])
            myHigh=int(sStuffH['binHighCount'])
            if myLow > myMid*BIN_DIFFERENCE_MID_MULTIPLIER:
                if myHigh > myMid*BIN_DIFFERENCE_MID_MULTIPLIER:
                    ## Ignore S.pal. where coverage > COVERAGE_SPAL_DEPTH_CUTOFF
                    #if s == "spal":
                    #    if sStuffH['avgCovDepth'] > COVERAGE_SPAL_DEPTH_CUTOFF:
                    #        if DEBUG: print "[binCounts] Ignoring coverage >",COVERAGE_SPAL_DEPTH_CUTOFF,":",g,s,sStuffH['avgCovDepth']
                    #        numSpalIgnoreH[g]=1
                    #        continue

                    testH[g][s]['inversedBins']+=1
                    numInversedBinsAllH[g]=1
                    
            #if g=="SPU_024901":
                #print "[binCounts]",g,s,testH[g]['numInversedBins'],myLow,myMid,myHigh
        
    print "[binCounts] Found",len(numInversedBinsAllH),"genes with inversed allele frequency bin counts low,high >",BIN_DIFFERENCE_MID_MULTIPLIER,"*middle bin."
    #print "[binCounts] (Ignored",len(numSpalIgnoreH),"S.pal. sequences where coverage >",COVERAGE_SPAL_DEPTH_CUTOFF,".)"
    
    return testH

#### Give some usage information for this script #######################################
def usage(errorNum):
    print "paralogTest.py - Score paralog"
    print 
    print "usage: paralogTest.py [-vDh] -a < alignment stats CSV >"
    print
    print "     Expected CSV columns: gene,species,avgCovDepth,biMorphCount,triMorphCount,binLowCount,binMidCount,binHighCount"
    print
    
    sys.exit(errorNum)


#### main #######################################
def main(argv):

    alignStatCsv=""
    
    try:
        opts, args = getopt.getopt(argv, "hD2:a:", ["help","debug","alignStatCsv"])
    except getopt.GetoptError:
        usage(20)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(21)
        if opt in ("-a", "--alignStatCsv"):
            alignStatCsv = arg
        if opt in ("-D", "--debug"):
            global DEBUG
            DEBUG=1
            
    if alignStatCsv == "":
        usage(33)

    geneH={}
    testH={}

    ## Expected format:
    #gene        spp     averageCoverage  biMorphCount  triMorphCount  binLowCount  binMidCount  binHighCount
    #SPU_009691  sdro    57.8281228669    326           0              46           556          46
    #SPU_009691  afraSS  56.00230571      246           0              57           373          56
    #SPU_009691  spal    16.6795851528    244           0              117          287          115
    #SPU_009691  sfra    51.1422319475    29            0              4            46           4
    #SPU_009691  hpul    63.8383346851    82            0              30           104          30
    #SPU_009691  snud    63.1647457627    25            0              13           24           13
    #SPU_009691  sprp    42.2957556096    87            0              31           112          31
    #SPU_009691  sintSS  58.1151662612    85            0              23           126          22
    #SPU_009691  pdep    79.904396658     14            0              4            22           5

    csvFP=open(alignStatCsv,"r")
    for line in csvFP:
        lineA=line.split(',')
        if lineA[0] == "gene":
            continue
        
        g=lineA[0]
        if not geneH.has_key(g):
            geneH[g]={}
            testH[g]={}
        
        s=lineA[1]
        geneH[g][s]={}
        
        geneH[g][s]['avgCovDepth']     =float(lineA[2])
        geneH[g][s]['biMorphCount']    =int(lineA[3])
        geneH[g][s]['triMorphCount']   =int(lineA[4])
        geneH[g][s]['binLowCount']     =int(lineA[5])
        geneH[g][s]['binMidCount']     =int(lineA[6])
        geneH[g][s]['binHighCount']    =int(lineA[7])
    
        #print "[main]",g,s,geneH[g].keys()
        
    csvFP.close()
    
    print "Found",len(geneH),"genes."
    testH,sppCovH = coverageScore(testH,geneH)
    testH = trimorphicCounts(testH,geneH)
    testH = binCounts(testH,geneH)

    numSpalIgnoreH={}
    testGoldH={}
    testSilverH={}
    
    silverCount=0
    outFP=open("paralog_test_results.summary.csv","w")
    outFP.write("gene,species,meanCovSppAllGene,stdCovSppAllGene,meanCovSppThisGene,coverageGt3SD,triMorphCount,inversedBins,binLowCount,binMidCount,binHighCount,goldTest\n")
    for g,mySpH in testH.iteritems():
        #print "[main]",g,mySpH
        all_isCov=0
        all_numTri_gt1=0
        all_numTri_gt2=0
        all_isInv=0
        silverScoresH={}
        silverScoresH["silver-"]=0
        silverScoresH["silver+"]=0
        silverScoresH["silver++"]=0
        for s,myH in mySpH.iteritems():
            lineA=[]
            lineA.append(g)
            lineA.append(s)
            #print "[main]",g,s,myH
            
            isCov=myH['coverageGt3SD']
            if DEBUG: print "[main] isCov",g,s,isCov
            if isCov >=1:
                all_isCov+=1

            numTri=myH['triMorphCount']
            if DEBUG: print "[main] numTri",g,s,numTri

            if numTri >=1:
                all_numTri_gt1+=1
            if numTri >=2:
                all_numTri_gt2+=1
                
            isInv=myH['inversedBins']
            if DEBUG: print "[main] isInv",g,s,isInv
            if isInv >=1:
                all_isInv+=1
            
            lineA.append(str(sppCovH[s]['covAvg']))
            lineA.append(str(sppCovH[s]['covStd']))
            lineA.append(str(geneH[g][s]['avgCovDepth']))
            lineA.append(str(isCov))
            lineA.append(str(numTri))
            lineA.append(str(isInv))
            lineA.append(str(geneH[g][s]['binLowCount']))
            lineA.append(str(geneH[g][s]['binMidCount']))
            lineA.append(str(geneH[g][s]['binHighCount']))
            
            ## Test Silver - just for this species
            silverScore="silver-"
            testSilverH[g]=silverScore
            
            ## "if a trimorphic call is present in Spal and coverage > 25, it is flagged as containing a paralog.
            ##  Otherwise, the trimorphic call is ignored." - GP
            if numTri>=1:
                if s == "spal":
                    if geneH[g][s]['avgCovDepth'] < COVERAGE_SPAL_DEPTH_CUTOFF:
                        if DEBUG: print "[main] ignoring S.pal. tri-morph due to low coverage",COVERAGE_SPAL_DEPTH_CUTOFF,g,s,geneH[g][s]['avgCovDepth']
                        numSpalIgnoreH[g]=1
                    else:
                        silverScore="silver+"
                        testSilverH[g]=silverScore
                else:
                    silverScore="silver+"
                    testSilverH[g]=silverScore
                
            if (isCov>=1) and (numTri>=1):
                silverScore="silver++"
                testSilverH[g]=silverScore

            if (isCov>=1) and (numTri>=1) and (isInv>=1):
                silverScore="silver++"
                testSilverH[g]=silverScore
            
            silverScoresH[silverScore]+=1
            lineA.append(silverScore)
                
            lineS=",".join(lineA)
            outFP.write(lineS+"\n")
        
        if ( silverScoresH["silver+"]>0 ) or ( silverScoresH["silver++"]>0 ):
            silverCount+=1
            if silverScoresH["silver+"]>0:
                testSilverH[g]="silver+"
            
            if silverScoresH["silver++"]>0:
                testSilverH[g]="silver++"
                
        ## Test Gold - any two species combination
        testGoldH[g]="gold-"
        
        ## Gold A. Alignments flagged as "gold+" would have TMC >=1 for >=2 species.
        if (all_numTri_gt1>=2):
            testGoldH[g]="gold+"
            
        ## Gold B. Alignments flagged as "gold+" would have TMC >=1 for =1 species
        ##      and inversed bin counts for >=2 species
        ##      (bin counts need not be elevated in the species with the TMC)
        if (all_numTri_gt1>=2) and (all_isInv>=2):
            testGoldH[g]="gold+"
            
        ## Gold C. Alignments flagged as "gold+" would have TMC >=1 for =1 species
        ##      and >2SD for >=2 species
        ##      (coverage not necessarily elevated in the species with the TMC).
        if (all_numTri_gt1>=2) and (all_isCov>=2):
            testGoldH[g]="gold+"
        
        ## Gold D. Alignments flagged as "gold++" would have TMC >=2 for >=2 species.
        if (all_numTri_gt2>=2):
            testGoldH[g]="gold++"
            
        ## Gold E. E. Alignments flagged as "gold+++" would have TMC >=2 for >=2 species
        ##      and inversed bin counts for >=2 species
        ##      and >2SD for >=2 species.
        if (all_numTri_gt2>=2) and (all_isInv>=2) and (all_isCov>=2):
            testGoldH[g]="gold+++"
            
        if DEBUG:
            print "[main]",g,"all_numTri_gt1=",all_numTri_gt1,"all_numTri_gt2=",all_numTri_gt2,"all_isCov",all_isCov,"all_isInv=",all_isInv,testSilverH[g],silverScoresH,testGoldH[g]
        
            
    outFP.close()
    print "Found", round(float(silverCount)/len(geneH),4)*100,"% of genes as likely paralogs per any 'Silver standard':",silverCount,"of", len(geneH)
    print "Ignored",len(numSpalIgnoreH),"S.pal. tri-morph due to low coverage <",COVERAGE_SPAL_DEPTH_CUTOFF
    
    scoreFP=open("paralog_test_results.scored.csv","w")
    scoreFP.write("gene,goldTestResult,silverTestResult\n")
    goldCount=0
    for g,score in testGoldH.iteritems():
        lineA=[]
        lineA.append(g)
        lineA.append(score)
        lineA.append(testSilverH[g])
        lineS=",".join(lineA)
        scoreFP.write(lineS+"\n")
        if score != "gold-":
            goldCount+=1
            
    scoreFP.close()
    
    print "Found", round(float(goldCount)/len(geneH),4)*100,"% of genes as likely paralogs per 'Gold standard':",goldCount,"of", len(geneH)
    
    outFP.close()
        

    
#### Start here. #######################################
if __name__ == "__main__":
    main(sys.argv[1:])
