# MARSS
Multiple Alignment of Reference and Short readS

MARSS is a tool developed to generate Multiple Species Consensus Alignments and quality control statistics for comparative genomics analyses of regions of a reference genome from aligned reads. It also implements a test to identify and score potential paralogs in whole genome or transcriptome sequencing comparative genomics studies. MARSS was developed with high throughput in mind and can be implemented in high performance computing environment for compatative genome analyses using next-generation short read technology.

marss - marss is a program to generate consensus coding sequences from whole-genome reference alignments of short read sequences.  

marssCodon - marssCodon is a program to generate consensus coding sequences from whole-genome reference alignments of short read sequences based on codons in phase.

## Requirements
     - Python 2.7
     - Python modules : pysam, BioPython, twobitreader, pp
     - Kent src tools

## Installation

Refer to INSTALL.txt at the project home page for help with installation of MARSS as well as required software.

## marss

### Usage

#### Input
- Input for marss is a comma-delimited list of the name and path to BAM file for each taxa (one per line), the path to the reference genome, and the path to the gene model. 
- marss expects the reference genome to be in 2bit format. 
- Reads must have previously been aligned to this reference and the alignments are provided as a BAM file (e.g., consistent coordinates and names of contigs in the genome 2bit and BAM files).
- The single gene model must be provided in the Gene Pred format.
- Program parameters are provided at runtime from command line switches. A description of all software options is available from the command line through the help output of the program.
- For more information on the format of the input files, reference README.DATA.txt

#### Options
     -2, --genome2bit         : reference genome 2bit file https://genome.ucsc.edu/FAQ/FAQformat.html#format7
     -p, --genPredFile        : gene prediction file (UCSC genPred table format) Https://genome.ucsc.edu/FAQ/FAQformat.html#format9
     -b, --bam ListFIle       : file containing comma delimited list of taxa name and path to the alignment BAM files. One taxa per line BAM indexes must be available as <file.bam>.idx http://samtools/github.io/hts-specs/SAMv1.pdf
     -m, --method             : method to use for consensus call [described below in Tabel 2 in Reference Table section]
     -a, --minorAllelCount    : minimum minor allele count required for an allele to be considered valid in the alignment and included (default 8)
     -f, --minorAllelFreq     : minimum minor allele frequency required for an allele to be considered valid in the alignment and included (default 0.126)
     -q, --minMapQuality      : minimum read mapping quality to include (default 25)
     -s, --minBaseQuality     : minimum base phread quality score to include (default 25)
     -X, --removeStopCodon    : truncate output sequences by 3 bases for codeml
     -P, --ignoreProperPair   : include all reads (default: ignore pairing status)
     -T, --test               : preform reporting of prerequisite packages
     -D, --debug              : set verbose output. You probably do not want to do this
     -h, --help               : shows this help (default: -h works with no other options)

#### Methods for nucleotide consensus calling in marss
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


#### Alignmnet statistics collected by marss for each taxa in the gene alignment
| Characteristic | Data Type | Description |
| :------------- | :-------- | :---------- |
| siteCount      | Integer   | The count of sites evaluated for the gene model |
| baseCalls      | Integer   | The count of successful base calls evaluated for the gene model (including no calls) |
| totalBases     | Integer   | The count of all of the bases evaluated |
| averageCoverage | Float    | baseCalls / count of nucleotides evaluated (i.e., gebe model length) |
| minorAlleleCalls | Integer | The count of alleles considered to be valid based upon frequency |
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
taxa    readName                                indelSite
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
>ref
ATGGATGAAGTCTTGAAACTTCTTCATAAACTCAATGTTAACCTCAGCAAGAGGAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>sdro
ATGGATGAAGTCTTGAAACTTCTTCAYAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>afraSS
ATGGATGAAGTCTTGAAACTTCTTCATAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>spal
ATGGAYGAAGTCTTGAARCTTCTTCAYAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAGCTCTTCAGGGAAGCCGATACCAACATCGATGAGCAC
>sfra
ATGGATGAAGTCTTGAAACTTCTTCACAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCYGATACCAACATCGATGAGCAC
>hpul
ATGGATGAAGTCTTGAAACTTCTTCACAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTCTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>snud
ATGGATGAAGTCTTGAAACTTCTTCACAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATTGATGAGCAC
>sprp
ATGGATGAAGTCTTGAAACTTCTTCATAAACTCAATGTTAACCTCAGCAAGAGRAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATMGAYGAGCAC
>sintSS
ATGGATGAAGTTTTGAAACTTCTTCACAAACTCAATGTCAACCTCAGCAAGAGAAAAGTCAAACAACTCTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>pdep
ATGGATGAAGTCTTGAAACTTCTTCACAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCAGATACCAACATTGATGAGCAC
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

## marssCodon

### Usage

#### Input
- Input for marssCodon is a comma-delimited list of the name and path to BAM file for each taxa (one per line), the path to the reference genome, and the path to the gene model. 
- marssCodon expects the reference genome to be in 2bit format. 
- Reads must have previously been aligned to this reference and the alignments are provided as a BAM file (e.g., consistent coordinates and names of contigs in the genome 2bit and BAM files).
- The single gene model must be provided in the Gene Pred format.
- Program parameters are provided at runtime from command line switches. A description of all software options is available from the command line through the help output of the program.
- For more information on the format of the input files, reference README.DATA.txt

#### Options
     -2, --genome2bit           reference genome 2bit file
                                https://genome.ucsc.edu/FAQ/FAQformat.html#format7
     -p, --genePredFile         gene prediction file (UCSC genePred table format)
                                https://genome.ucsc.edu/FAQ/FAQformat.html#format9
     -b, --bamListFile          file containing comma delimited list of taxa name and
                                path to the alignment BAM files. One taxa per line.
                                BAM indexes must be available as <file.bam>.idx
                                http://samtools.github.io/hts-specs/SAMv1.pdf
     -m, --method               heterozygote site method to use for consensus call [described below]
     -n                         number of CPUs/cores to use for parallel BAM file processing (1)
     -g                         evaluate for glycan motifs (N*S or N*T) (experimental and memory intenstive)
     -D, --debug                set verbose output. You probably do not want to do this.
    (-h) --help                 show this help (-h works with no other options)

#### Methods for nucleotide consensus calling in marss
| Method | Description |
| :----- | :---------- |
|freq|majority observed codon|
|ref|reference codon|

## Output

#### Description of output files
| Output File Name | Description |
| :--------------- | :---------- |
| alignInfo.txt    | Reports various summary information for the alignment per taxa |
| consensAlign.fa  | Reports the MSCA |
| obsCodons_\<taxa\>.txt | Reports the called and counts of each codon at each position for "taxa".|
| glycanMotif_\<taxa\>.txt | Reports the called and counts of each motif at each position for each read of "taxa".|

#### Alignment statistics collected by marssCodon for each taxa in the gene alignment
| Characteristic | Data Type | Description |
| :------------- | :-------- | :---------- |
| codonCount     | Integer   | The count of codon sites evalauted for the gene model |
| hetSitesConserv   | Integer   | The count of heterozygote sites identified that match the reference |
| hetSitesReplace   | Integer   | The count of heterozygote sites identified that differ from the reference |
| hetSitesSilent   | Integer   | The count of heterozygote sites identified that are silent replacement (i.e., nonsynonymous) |
| nnnCount       | Integer   | The count of sites identified that were uncalled |
| nnnSites       | Integer   | The position of successful base calls evaluated for the gene model (including no calls) |
| polyAlleleCount   | Integer   | The count of polymorphic sites |
| polyAlleleSites   | Integer   | The positions of the polymorphic sites |
| chosenCodon | String | The codon selected by the method choosen at the position |
| obsCodonCounts | String | The counts for each codon observed at the position |
| motifCount_NXS   | Integer   | The count of sites having the glycan motif (N\*S)  |
| motifCount_NXT   | Integer   | The count of sites having the glycan motif (N\*T)  |
| motifSites_NXS   | Integer   | The positions of sites having the glycan motif (N\*S)|
| motifSites_NXT   | Integer   | The positions of sites having the glycan motif (N\*T) |
| motifSites_glycanSiteCount   | Integer   | The count of glycan motifs at at each position |
| motifSites_glycanSiteFreq   | Float   | The frequency of sites having a glycan motif |
| motifSites_glycanSiteMask   | String   | A mask indicating the presence or absense of a detected glycan motif at the position |
| glycanMotif   | String   | The predicted glycan motif |
| glycanMotifDetail   | String   | Detailed glycan motif |
| motifScore | Float | Score for the motif call |

#### Examples of Outputs

obsCodons_mytaxa.txt:
<pre>
179     12816   2       TGC     TGC:41
180     12819   2       AGA     AGA:38
181     12822   2       TGC     TGC:22,TGT:17
182     12825   2       GTG     GTG:38
183     12828   2       GAA     GAA:40
184     13819   3       NNN
185     13822   3       GAC     GAC:26
186     13825   3       ACA     ACA:27
187     13828   3       TGG     TGG:27
188     13831   3       GAT     GAT:27
</pre>

consensAlign.fa :
<pre>
>ref
ATGGATGAAGTCTTGAAACTTCTTCATAAACTCAATGTTAACCTCAGCAAGAGGAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>sdro
ATGGATGAAGTCTTGAAACTTCTTCAYAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>afraSS
ATGGATGAAGTCTTGAAACTTCTTCATAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>spal
ATGGAYGAAGTCTTGAARCTTCTTCAYAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAGCTCTTCAGGGAAGCCGATACCAACATCGATGAGCAC
>sfra
ATGGATGAAGTCTTGAAACTTCTTCACAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCYGATACCAACATCGATGAGCAC
>hpul
ATGGATGAAGTCTTGAAACTTCTTCACAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTCTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>snud
ATGGATGAAGTCTTGAAACTTCTTCACAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATTGATGAGCAC
>sprp
ATGGATGAAGTCTTGAAACTTCTTCATAAACTCAATGTTAACCTCAGCAAGAGRAAAGTCAAACAACTGTTTAGGGAAGCCGATACCAACATMGAYGAGCAC
>sintSS
ATGGATGAAGTTTTGAAACTTCTTCACAAACTCAATGTCAACCTCAGCAAGAGAAAAGTCAAACAACTCTTTAGGGAAGCCGATACCAACATCGATGAGCAC
>pdep
ATGGATGAAGTCTTGAAACTTCTTCACAAACTCAATGTTAACCTCAGCAAGAGAAAAGTCAAACAACTGTTTAGGGAAGCAGATACCAACATTGATGAGCAC
</pre>

glycanMotif_mytaxa.txt:
<pre>
codonStart     codonEnd  exonNumber     siteStart      readName                                glycanMotif    glycanMotifDetail  motifScore
14             16        0              10047          HS2:148:C0EN2ACXX:4:1205:10328:24507    NLS            N14L15S16           1
14             16        0              10047          HS2:148:C0EN2ACXX:4:1107:14271:75957    NLS            N14L15S16           1
14             16        0              10047          HS2:148:C0EN2ACXX:4:2101:7979:93625     NLS            N14L15S16           1
14             16        0              10047          HS2:148:C0EN2ACXX:4:1307:3613:59334     NLS            N14L15S16           1
14             16        0              10047          HS2:148:C0EN2ACXX:4:2102:16618:89925    NLS            N14L15S16           1
14             16        0              10047          HS2:148:C0EN2ACXX:4:2102:3860:155303    NLS            N14L15S16           1
14             16        0              10047          HS2:148:C0EN2ACXX:4:2301:18438:183010   NLS            N14L15S16           1
14             16        0              10047          HS2:148:C0EN2ACXX:4:1105:16887:18933    NLS            N14L15S16           1
14             16        0              10047          HS2:148:C0EN2ACXX:4:2208:2613:54804     NLS            N14L15S16           1
</pre>

## Testing
Refer to the DEMO.txt file at the project home page (http://github.net/kordk/marss) in order to see how to run a demo of marss and marssCodon.

## Examples of studies which used marss and marssCodon

Kober KM, Bernardi G. Phylogenomics of strongylocentrotid sea urchins. BMC Evolutionary Biology. 2013. 13:88. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3637829"> PMC3637829</a>.

Kober KM, Bernardi G. Erratum to: Phylogenomics of strongylocentrotid sea urchins. BMC Evol Biol. 2017 Feb 13; 17(1):50. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5307700">PMC5307700</a>.

Kober KM, Pogson GH. Genome-wide signals of positive selection in strongylocentrotid sea urchins. BMC Genomics. 2017 07 21; 18(1):555. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5521101"> PMC5521101</a>.

Saarman NP, Kober KM, Simison WB, Pogson GH. Sequence-Based Analysis of Thermal Adaptation and Protein Energy Landscapes in an Invasive Blue Mussel (Mytilus galloprovincialis). Genome Biol Evol. 2017 Oct 01; 9(10):2739-2751. <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5647807"> PMC5647807</a>.

## Authors
Kord Kober : kord.kober@ucsf.edu

Samantha Danison : samanthadanison00@gmail.com

Grant Pogson :  pogson@ucsf.edu

## Acknowledgements

Funding for sea urchin data collection was provided to Dr. Kober by the National Science Foundation (DEB-1011061), the STEPS Foundation, Friends of Long Marine Lab, and the Myer’s Trust. The funding bodies did not participate in the design of the study or collection, analysis, and interpretation of data or in writing the software. We are extremely grateful to Giacomo Bernardi (UCSC) for advice and guidance on the development and evaluation of the evolutionary analyses and Satish Pillai (UCSF) for advice and guidance on the development of the glycan testing module in marssCodon.
