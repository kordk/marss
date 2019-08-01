# MARSS
Multiple Alignment of Reference and Short readS

MARSS is a tool developed to generate Multiple Species Consensus Alignments and quality control statistics for comparative genomics analyses of regions of a reference genome from aligned reads. MARSS was developed with high throughput in mind and can be implemented in high performance computing environment for compatative genome analyses using next-generation short read technology.

## Requirements
     - Python 2.7
     - Python modules : pysam, BioPython, twobitreader, Kent src tools

## Installation

Refer to INSTALL.txt at the project home page for help with installation of MARSS as well as required software

## Usage

- Input for MARSS is a comma-delimited list of the name and path to BAM file for each taxa (one per line), the path to the reference genome, and the path to the gene model. 
- MARSS expects the reference genome to be in 2bit format. 
- Reads must have previously been aligned to this reference and the alignments are provided as a BAM file (e.g., consistent coordinates and names of contigs in the genome 2bit and BAM files).
- The single gene model must be provided in the Gene Pred format.
- Program parameters are provided at runtime from command line switches. A description of all software options is available from the command line through the help output of the program.

#### Options
     -2, --genome2bit : reference genome 2bit file https://genome.ucsc.edu/FAQ/FAQformat.html#format7
     -p, --genPredFile : gene prediction file (UCSC genPred table format) Https://genome.ucsc.edu/FAQ/FAQformat.html#format9
     -b, --bam ListFIle : file containing comma delimited list of taxa name and path to the alignment BAM files. One taxa per line BAM indexes must be available as <file.bam>.idx http://samtools/github.io/hts-specs/SAMv1.pdf
     -m, --method : method to use for consensus call [described below in Tabel 2 in Reference Table section]
     -a, --minorAllelCount : minimum minor allele count required for an allele to be considered valid in the alignment and included (default 8)
     -f, --minorAllelFreq : minimum minor allele frequency required for an allele to be considered valid in the alignment and included (default 0.126)
     -q, --minMapQuality : minimum read mapping quality to include (default 25)
     -s, --minBaseQuality : minimum base phread quality score to include (default 25)
     -X, --removeStopCodon : truncate output sequences by 3 bases for codeml
     -P, --ignoreProperPair : include all reads (default: ignore pairing status)
     -T, --test : preform reporting of prerequisite packages
     -D, --debug : set verbose output. You probably do not want to do this
     -h, --help : shows this help (default: -h works with no other options)

#### Methods for nucleotide consensus calling in MARSS
| Method | Description |
| :----- | :---------- |
| mda | Polymorphic sites are called as the majority divergent allele as compared to the reference |
| ma | Polymorphic sites are called as the majority allele observed from valid reads |
| poly | Polymorphic sites are called as polymorphic using the IUPAC nomenclature |
| ref | Polymorphic sites are called the reference if observed, else called the majority allele |
| hetRefAll | Polymorphic sites are called reference regardless if the reference allele is observed |
| fixed | Polymorpic sites are randomly called from the observed nucleotides |


## Output

#### Description of output files
| Output File Name | Description |
| :--------------- | :---------- |
| indelInfo.txt    | Reports the sites of insertion or deletion for each taxa and readName |
| alignInfo.txt    | Reports various summary information for the alignment per taxa |
| consensAlign.fa  | Reports the MSCA |
| siteInfo.txt | Reports the frequency of each allele and the siteCoverageDepth at each reference position of each taxa |


#### Alignmnet statistics collected by MARSS for each taxa in the gene alignment
| Characteristic | Data Type | Description |
| :------------- | :-------- | :---------- |
| siteCount      | Integer   | The count of sites evaluated for the gene model |
| baseCalls      | Integer   | The count of successful base calls evaluated for the gene model (including no calls) |
| totalBases     | Integer   | The count of all of the bases evaluated |
| averageCoverage | Float    | baseCalls / count of nucleotides evaluated (i.e., gebe model length) |
| minorAlleleCalls | Integer | The count of alleles considered to be valid based apon frequency |
| divergentAlleleCalls | Integer | The count of alleles different from the reference |
| polySiteCount | Integer | The count of valid polymorphic sites observed |
| biMorphCount | Integer | The count of polymorphic sites wirh 2 valid alleles observed |
| triMorphCount | Integer | The count of polymorphic sites with >2 valid alleles observed |
| delCount | Integer | The count of deletions observed |
| indelCount | Integer | The count of insertions observed |
| multimapCount | Integer | The count of read multi-mappings observed |
| indelSite | Integer | The site of intersection |

#### Examples of Outputs


indelInfo.txt :
<pre>
taxa   readName                                 indelSite
sdro    HS3:171:d0le4acxx:1:1105:18332:137157   28503 
sdro	HS3:171:d0le4acxx:1:1205:16883:73266	32073 
sdro	HS3:171:d0le4acxx:1:1108:16497:10621	28503 
sdro	HS3:171:d0le4acxx:1:1307:9119:139026	28503 
sdro	HS3:171:d0le4acxx:1:1307:18934:36078	28503 
sdro	HS3:171:d0le4acxx:1:2203:15217:27563	34864 
sdro	HS3:171:d0le4acxx:1:1206:10535:154891	32069 
sdro	HS3:171:d0le4acxx:1:1206:10535:154891	32073 
sdro	HS3:171:d0le4acxx:1:1106:6061:52179	28503 
sdro	HS3:171:d0le4acxx:1:1102:3279:38879	28503
</pre>

alignInfo.txt :
<pre>
taxa    baseCalls   delCount   polySiteCount   multimapCount   totalBases   divergentAlleleCalls   averageCoverage   triMorphCount   siteCount   biMorphCount   minorAlleleCalls   indelCount 
sdro	7325        0          326             828             378878       7651                   57.8281228669     0               7325        326            7651               113 
afraSS	7373        0          246             1136            360892       7618                   56.00230571       0               7373        246            7618               200 
spal	7328        0          244             284             117733       7599                   16.6795851528     0               7328        244            7599               42 
sfra	7312        0          29              136             184669       7328                   51.1422319475     0               7312        29             7328               40 
hpul	7398        0          82              1487            424232       7479                   63.8383346851     0               7398        82             7479               149 
snud	7375        0          25              1617            395043       7384                   63.1647457627     0               7375        25             7384               191 
sprp	7398        0          87              3002            276171       7485                   42.2957556096     0               7398        87             7485               47 
sintSS	7398        0          85              1580            373686       7480                   58.1151662612     0               7398        85             7480               181 
pdep	7301        0          14              2167            512438       7312                   79.904396658      0               7301        14             7312               261
</pre>

consensAlign.fa :
<pre>
sdro
ATGGATGAAGTCTTGAAACTTCTTCAYAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATCGATGAGCACTTGGGCAAGCTGGACTTTGATGAGTTTGTGCACTTCTACAAGACATTGTCGATGCGCCCTGAGCTGTATGGGCTTCTCAGAGAGTACAGCTGTGGCAAGGAACACATGACGCTGGAGGACGTGGATGTGTTTCTAAGGAATGAACAGGATATTAACAAACTCAATGAGGAGGGATGTAATAGCCTCATTGAAAAGTATGAACCGGTACCTGAGAATATCTGTACTGGCAGGCTGGGCATTGATGGACTGATGAGGTACCTCCTGAGTGAGGAGGGTGACCTTTTTAACTCAGGTCATCGAGGGGTCAACCAGGACATGACCCAGCCCCTTTCACACTACTACATCGCCTCTTCACACAACACGTACTTACTCGGTGACCAGCTGATGTCACAATCTAGTGTAGATGTCTATGGTCTCGTCCTGCAGGCGGGCTGCAGATGTGTGGAATTGGACACATGGGATGGGAAGGACGGTGAACCGGTCATCTACCATGGTTACACGCTGACCACAAAGATCAAGTTCAGAGATGTCATCACAGTGGTCAACAAGTATGCCTTTATGTCCAGTCCCTACCCTGTCATCCTGTCCATAGAGAATCACTGTACCTTAGAGCAGCAGAAGAAGATGGCTAGATACCTTCTAGAGATCCTTGGTGATCAACTCATCATATCAAGCCCACCTGAAGGGAGTTCAGGAGGACTTCCTTCTCCAGAATCACTCAAGTATAAAATACTTATCAAGGCMAAGAAGCTTCCATTAGACCATGATGCTGATCTCGAGGCTGGTGAAGTTAGTGAAGAAGATAGTGCTGATGAACTTGATGAAGACTTCAAGTTGGAGAAGACCGAGACGAGAAGTAAGTTTGAATCTATTGCCATGGCACAGCTTGCATTGATGAGGAAGAAACCTATCTCACCATCGCAACAGGCATGGCAGAAGATAGCGCACAAGCTCAAAGAAAAGAACACTGAAAATAAAGCTTCCAARAAYGGGAAGACTGGACTCACAGGAACCTTCAGGAAATTAAAAGCTACMAGAAGGMGWGCCAAGCTGCAGCACTCGTCCGAGTCAGATAACAACAGCAGTATGGAATACGAKCGATCATCTTCTCGAGATGATGACTTAGAGGGCGGTGCCAATGAGAACCAGCAGGCAGAGCAAGCGAAACAACGGATGAGCTCAGCYAACAGAAAGAGGGCGTTTGTTCTCTCAAGACAGYTGTCAGATCTAGTCAAGTACACCAGATCGGTTGCATTCAAGGGCTTCCCAGAAAACATGCCATGTTGGGAGCTGCCCTCGCTTGGTGAGATGAGRGCAAACAATCTTGCGATGACGAGAGCTGCAGACTTTGCTACCTTCAGTACTCACAACATGTGCAGAGTYTACCCATCGGCTTATCGGATAGACTCCAGTAACTTCAATCCTCAGCCGCTATGGAACTGTGGATGCCAGCTTGTTGCTCTGAACTACCAGACTGAAGGTCGTGCTATGCAATTGTTGAGAGCTAAGTTYAGAGCCAATGGGAGYTGTGGTTATGTGCTCAAACCTAGCATCCTCAGACAAGACAAAGTCCATTTTGATGCAACATCCAATAATTCCATTGCTGGCGTCCAGAAGAAGCAACTGACYATTCACATYATCAGCGGTCARCAGCTACCCAAACCACCMCAAAGCCTGCTGGGAGAACGTGGAGAGATCATTGACCCCTATGTRGAAGTKGAAGTWGTTGGTCTGAATGGWGATTGCACCAAGGCACAGACTAGGACTGTYCAAGATAAYGGTTTCAATCCAGTCTGGGACTATGTWGTGACGTTCCCCGTCACTCTGCCYGAGCTTACYCTGGTCAGGTTTGTWGTATGGGATGAAGATCCCATCGGKAGAGAYTTCATTGGCCAGGCAACCTTYCCTTTCACCAGCCTTTGTCAAGGGTACCGTCACATCCATCTTGAGGGTCTTGACCATGCAACTGTCTTTGTTCATGTCACCATTGAAGACTTCACTGATAAGAAATACATCAAAAACAGCAAAGGCCTGCGCTCCTCCTCCATCACCCGCACCACTTCCTCCCCGAGCCGCAACCGAACCAAATCCGTTGACACGGAGGTCAATCTGATGACCGGTGGAAGCCCCTCCCATCAGGCCTCGCCCAAGAGACGCACCTCCCTGGCAGAGAGGCTTGGTATCAGGAGGAGACACAGCACCACAGCACTCCTGGTGGAAGGCATCATGAATGTGAGYRGTGGCASCGGTGGGATCAAGAGTCCTACCTCACCCACTTCCACACCCACSTTTGGTCTGCTGAGGCAGAAGAGAAGGARCTATGAAGGAGAGATGAGGAAGGGTCGCTCACAAGGTGATGAGGAAGGATGCGATGATGTTTTTGGKGATGACGACTGTGACTTGAGTCCTGTTGACTGTGATGCCATTTTGGATCGCATCATGTCCAGTGAGGTAGCSGTAGAGCTAGAGGCWGATGGTAGTGTGAGTGTTCCTCTCAATGAGATTATWCACATCCTCCGCCTAGGAATTTCTCCTTCACAAATCATCCAAGTCATCCAAGATGGCCGCTTGAGTCTGAGACCAGCRCTCCATGAGGTTCCTCCCATTAGYTCAGTAGGAKAKTCTTTGGATGTTGCTGCTARARGTGAACTAGAGTCCAGTCYAATGGAYACTCAAGAAGGATTAGGTAGYGGAATCCTTTCCAAGTTAGAACATTWTGACCMTGAGCCAGCTCCTACCTTGAGTARAATCTTAGACGATGATGATCTTGTGTTTCCTATGAATAATAATGATATCTCAACAGAGATGGAGGCCCTYCTGTCTGCYCGCTTYACGAAGGAGAGACTAGCCTTTGGGATGGTGTCSTATGATGAGTGTTCTTCGACTGAAGACTGCTCAAGTCTGGAATCTGCTGAGGTTTWTGTGACTGACACTCCAATAACKGACCnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnAACTGACCATGTGACCAATGACACATCAACTGACCATGTGACCGATAACAAATCAACTGAYCATGYGACCAATGACACATCAACTGACCATGYGACYAATRAYACATCAACTGAYCATGYGACCRATRACAMATCAACTGAYCATGTGACCAATRAYAMATCAASTGAYCATGTGACCRAWAATGCAAMAGGTGAACTCCTTTCAGRAGAGAATSGTCTTGTCTCTGGAGAGGTTCATGCAAAGCCAGAGCAAGAGATAGAAATACAAGAWAGAGAAATAGGAKTACCACARYTAWCAAACTCCTCTGAATACAAKGTAGAAGRTAATAGTRAAGCAGATGGGATTCTAGACGACCCCTGTCTCCAKGAAWATCAGGGTAGTGACTATGATAACCAGGRGGATATTGTGTATGATTTTGTAGAACAAGATGATAWTGTAGTACATGTACCTGGACAKGKTMARSTSGCTTCACAGCCAAAYGAATTTGTCCAAAGGTYRAATSAAGATGAYGTTTYTTYRAACAGTARTGAGGGTATGTGTCTGTCRAGCMWGGGGAATCCTARTAAATTMTCTTCATGTARTGGRAATAATGACGTAGTCRATGCATACAGTAATRCYWCCTATATMAGAAAGAATGARGATCATTATAATCTTAGCTCAGTAGATCACTCTCAAGKGTATATGCAATATCCAGAGGACAGAGAAGGGCGGTTATCWCAAGAGTGTGTCAAACAACTTGTACAAGGTCAGACAGAAGTAGAAMTRGATTCRGT
</pre>

siteInfo.txt :
<pre>
taxa	position   A               C                T      G     siteCoverageDepth
sdro	35236	   1.0	           0.0              0.0    0.0   67
sdro	35237	   0.0	           0.0              0.0    1.0   69
sdro	35234	   0.0	           1.0              0.0    0.0   66
sdro	35235	   0.608695652174  0.391304347826   0.0    0.0   69
sdro	35232	   0.0	           0.0              1.0    0.0   66
sdro	35233	   0.0	           1.0              0.0    0.0   68
sdro	35230	   0.0	           1.0              0.0    0.0	 63
sdro	35231	   0.0	           0.0              1.0    0.0   68
sdro	35548	   0.0	           1.0              0.0    0.0   42
</pre>

## Testing
Refer to the DEMO.txt file at the project home page (http://github.net/kordk/marss) in order to see how to run a demo of MARSS.

## Authors
Kord Kober :
kord.kober@ucsf.edu

Samantha Danison :
samanthadanison00@gmail.com

## Acknowledgements
TBD
