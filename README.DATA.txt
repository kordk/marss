The software tool expects three input files. Each file type is described
briefly below and includes an example for testing. That said, one file is a
list pointing to other files (i.e., BAM) used in the analysis.

1. Reference genome file. This will need to be formatted as a 2bit file.
    Usually, reference genomes are used as FASTA formatted files. In
    this case, the FASTA file will need to be converted into a 2bit
    file. This can be done using faToTwoBit from the Kent src tools (from
    the UCSC Genome Browser).

    E.g., To convert the reference FASTA strPur4.fa:

        faToTwoBit strPur4.fa strPur4.2bit

    Sample:
        http://hgdownload-test.cse.ucsc.edu/goldenPath/strPur4/bigZips/strPur4.2bit

3. Gene model. This will need to be formatted as a GenePred data file with an extra column.
    This is expected to be tab delimited export of a genePred table.
    (SEE: https://genome.ucsc.edu/FAQ/FAQformat.html#format9)
    The additional column is first, is labeled 'bin', and must contain
    an integer. It is a byproduct of the UCSC Table Browser export and
    is ignored by the tool, but expected as part of the format.
    Only one model is expected. 
    No header is expected.
    
    E.g., With headers for clarity:
bin name       chrom     strand txStart txEnd  cdsStart cdsEnd exonCount exonStarts            exonEnds
73  SPU_010036 Scaffold1 +      127266  237263 127266   237263 3         127266,144788,237142, 127296,144998,237263,

    E.g., Without headers as expected by the tool:
73  SPU_010036 Scaffold1 +      127266  237263 127266   237263 3         127266,144788,237142, 127296,144998,237263,

    Sample: This file will contain all of the gene models for this
        genome. Only one is used as input for the tool at a time.
    http://genome-preview.ucsc.edu/cgi-bin/hgTables?db=strPur4&hgta_group=genes&hgta_track=GLEAN_3_1_chado_UTR_gene&hgta_table=GLEAN_3_1_chado_UTR_gene&hgta_doSchema=describe+table+schema
                
3. List of alignment files. This comma-delimited file contains information for the
    alignments, one taxa per line. 

    E.g.,
    taxa1,/path/to/taxa1.bam
    taxa2,/path/to/taxa2.bam

4. Alignment files (one per taxa). This will need to be formatted as a
    BAM file with an index.

    E.g., Nine species were used in Kober & Pogson (2017). The alignment
        files are prohibitively large for sharing (189G total, 15-27G each).

    22G /mnt/frink/hilo/lantern/phd/proj/strPur4/sam/ssaha2.sdro_d206_KMK14.2013-03-19/mapped.merge.sorted.bam
    11G /mnt/frink/hilo/lantern/phd/proj/strPur4/sam/ssaha2.spal_p20a_KMK002.2013-04-21/mapped.merge.sorted.bam
    23G /mnt/frink/hilo/lantern/phd/proj/strPur4/sam/ssaha2.afra_af10_KMK12.2013-04-05/mapped.merge.sorted.bam
    27G /mnt/frink/hilo/lantern/phd/proj/strPur4/sam/ssaha2.sint_int1_KMK013.2013-04-25/mapped.merge.sorted.bam
    15G /mnt/frink/hilo/lantern/phd/proj/strPur4/sam/ssaha2.sfran_sf4_KMK010.2013-04-29/mapped.merge.sorted.bam
    24G /mnt/frink/hilo/lantern/phd/proj/strPur4/sam/ssaha2.snud_nu1_KMK011.2013-05-04/mapped.merge.sorted.bam
    22G /mnt/frink/hilo/lantern/phd/proj/strPur4/sam/ssaha2.pdep_pd3_KMK015.2013-04-17/mapped.merge.sorted.bam
    24G /mnt/frink/hilo/lantern/phd/proj/strPur4/sam/ssaha2.hpul_hp5_KMK16.2013-04-11/mapped.merge.sorted.bam
    21G /mnt/frink/hilo/lantern/phd/proj/strPur4/sam/ssaha2.spurp_SRA051151.2013-05-12/mapped.merge.sorted.bam

    Note that the BAM index should be in the same path as the BAM file.  E.g.,

        /path/to/taxa1.bam
        /path/to/taxa1.bam.bai

