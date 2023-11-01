# MacroHapCaller

* __phased multiple alleles calling for Next Generation sequencing generated long sequencing reads (designed for), e.g. PacBio HiFi__

*  __call variants of STR/SSR, SNP, InDel at the same time__

* __run faster with parallel computing__
  
* __also works for short reads__

Macrohaplotype caller for STR, SNP, and InDel from NGS long-read sequencing data, especially for PacBio HiFi reads
The MacroHapCaller calls targeted STRs/SSRs, SNPs, and InDels sites simultaneously in a row from each NGS read and clusters the variants into phased haplotype strings. MacroHapCaller is the best tool to analyze all haplotypic variants from genetically inherited DNA. It suits diploid, polyploid, and complex DNA mixtures from many individuals, e.g. DNA forensics. MacroHapCaller is programmed in Java with parallel computing enabled so it can run on any computing platform. 

The MacroHapCaller is also integrated into a pipeline with a shell script or Python3. The pipeline takes the BAM or fastq reads as input, align reads to the genome, sort, and index alignment, and then user can use the output to run MacroHapCaller. Additional tools are also available as util tools. For util usage : `java -jar appName.jar` where appName should be replaced with a detailed tool name. 


MacroHapCaller runs fast. It takes ~ 2 mins for merged data from four full runs of PacBio SMRTcell HiFi reads on a Linux machine with 12 threads, and a maximum allowed 120G RAM. 

## The Latest version
V0.4

## Historical versions
V0.3

## Environment
The Java run environment is needed, which is installed in most computer. You may just need to update it to the latest version. The Latest LTS Java 21 or higher is recommended. 

Step 1. Go to the Orcale website, download and unzip the file https://download.oracle.com/java/21/latest/jdk-21_linux-x64_bin.tar.gz ; 

Step 2. Add the path to your unzipped Java bin directory before running the command. e.g.
the Java is unzipped into c/java21, then the path will be c/java21/jdk-21/bin, the running command will be :

`path=c/java21/jdk-21/bin`

`$path/java -jar MacroHapCaller##.jar`

Or you can install the Java  JDK 21 in your computer, then you can run without $path/ before java , just like instructed below.

Or add the path into to environment.

## Installation:
The MacroHapCaller is programmed in Java and compiled. So just download the APP and run it.

### Install  MacroHapCaller from the Github. 
`git clone https://github.com/XuewenWangUGA/MacroHapCaller`

`cd MacroHapCaller`

Or just download all folders from Github and uncomppress the zip file. then go to the unzipped folder.

`cd MacroHapCaller`



### Install dependency (optional): fastq-filter 

 fastq-filter is the only needed dependency if want to filter the reads. Details on https://github.com/LUMC/fastq-filter. To install

`pip install fastq-filter`


## Help
Open a terminal window. Type the following command for help: 

`java -jar MacroHapCaller##.jar`

where ## is the version number, e.g. :   0.4.

##  Command and options: 
to change the version number as needed.

usage: `java -jar -Xmx100G MacroHapCaller0.4.jar [options]`

 **-a,--strAnchorFile** <arg>  &nbsp; configure file of anchors for STRs in
                                 tabular plain text
 
**-c,--MinReadCount** &nbsp; integer, minimum read count to report a
                                 macrohaplotype, default [2]
 
 **-d,--inDelBedFile** <arg> &nbsp; configure file of InDels in BED format in
                                 tabular plain text
 
 **-i,--input** <arg> &nbsp; input BAM file with a path
 
 **-l,--homopolymerLen** &nbsp; integer, homopolymer length threshold to
                                 skip call variant,default [4]
 
 **-m,--MaxMismatches**  &nbsp; integer, maximum distance of edit
                                 distance in anchor match, default [1]
 
 **-n,--snpPanelPosFile** <arg> &nbsp; configure file of SNPs in tabular plain
                                 text
 
 **-o,--output** <arg> &nbsp; output file name for macrohaplotype
                                 result
 
 **-p,--MinReadProportion** &nbsp; double, minimum read proportion of all
                                 reads for each locus to report a macrohaplotype, default
                                 [0.001]
 
 **-q,--LowQualThreshold** &nbsp; integer, phred-scale quality threshold
                                 for variant bases, default [15]
 
 **-r,--refGenomeFastaFile** <arg> &nbsp; Genome reference sequence file in .fasta
                                 format
 
 **-t,--ThreadNumber** &nbsp; integer, the number of computing threads,
                                 default [12]

 
 ## Example : 
 Download the subfolders and put them inside the folder of "MacroHapCaller". Then users can use parameters and configure files as used in the demo listed below.
 
    java -jar -Xmx120G MacroHapCaller0.4.jar -i demo_data/hg002.8kampl.Q30.4kreads.fastq.gz_GRCh38.bam -o hg002.8kampl.Q30.fastq.gz_GRCh38.HapVar.tsv -a config_v0.3/CODISSTR_anchor.XW.config_v0.3.txt -d config_v0.3/MHindels_v0.3.bed -n config_v0.3/MHsnps.pos_v0.3.txt -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta -l 2 -m 1 -q 15 -c 100 -p 0.01 -t 12

    the option -Xmx120G is to set allowed maximum memory allowed for this app. it can be ignored or changed to a small memory , e.g., -Xmx4G

## Build a pipeline from fastq reads or PacBio bam to macrohaplotypes


### Step 0: Prepare a genome reference file which should be the one used for alignment
E.g., human genome hg38: 
Download the genome sequence of human from the 1000 Genome Project to the folder "VarSeqStitcher" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
`wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa`

rename the genome sequence file as a short name:  

`mv GRCh38_full_analysis_set_plus_decoy_hla.fa genome.fasta`

   Then index the genome sequence with samtools (tool link: https://www.htslib.org/)
   samtools faidx genome.fasta
   
### Step 1: Preprocess pipeline for MacroHapCaller
 This pipeline will preprocess the original PacBio HiFi reads in BAM to fastq, read quality control, umi-analysis (optional), statistical summary, read-reference alignment, sort, and index. The major results are a read-reference alignment in the BAM format and BAM index, as well as statistical report files. Edit the path and files name to your specified names.
 
 `MH_umi_dedup_map_pipe.sh`

 Assume the alignment output will be named as "alignment.bam"
 

 ### Step 2: Call macrohaplotypes

 `java -jar -Xmx100G MacroHapCaller0.4.jar -i alignment.bam -o out.tsv -a config_v0.3/CODISSTR_anchor.XW.config_v0.3.txt -d config_v0.3/MHindels_v0.3.bed -n config_v0.3/MHsnps.pos_v0.3.txt -r genome.fasta -l 2 -m 1 -q 15 -c 100 -p 0.01 -t 12`
 
 Replace the real file name with your file name, and add file location path before the file if necessary. The above command uses the configure files coming with the software for human CODIS STR based HG38.
 
 The output will be in file: out.tsv, which can be viewed in spreadsheet, Mircosoft Excel or any other text file editor.
 
 
 

## Input of MacroHapCaller
 1. a read-reference alignment file in [BAM format](https://en.wikipedia.org/wiki/SAM_(file_format)). This alignment file could be generated from  `MH_umi_dedup_map_pipe.sh` or user-generated alignment file in .Bam format

 2. Configure files of targeted variant positional information for each of STRs, SNPs, and InDels

## Configure format for STR
The configure file listed the STR loci information in a tabular separated text file. An example file for forensic CODIS STR of human genome is given in "CODISSTR_anchor.XW.config_v0.2.txt"
Each line for one STR locus. There is one fixed headline before starting the detailed STR locus. Every line has 11 columns. The integer value for the last column is the estimated maximum length of TR, where 500 is enough for most cases. The coordinate starts from 1.

    Chrom	ChromStartPos_Str	ChromEndPos_Str	Name	Repeat_unit_length	Motif_must_present	Inner_offset	Anchor_left	Anchor_right	LeftAnchorStartPos	RightAnchorEndPos	MaxRefLength
    chr1	230765214	230765259	D1S1656@STR_118267	2	CACACA	0	CAGAAAATGAGAACACATG	GGTTATGCCAAAAGGGC	230765190	230765281	500
    chr1	230769616	230769683	D1S1656	4	TCTATCTATCTA	0	TCAGAGAAATAGAATCACTA	TGAGCAACACAGGCTTGAAC	230769556	230769721	500
    chr2	1489653	1489684	TPOX	4	AATGAATG	0	CAGAACAGGCACTTAGGGAA	AACGCTGACAAGGACAGAAG	1489624	1489716	500

    

## Configure format for InDel
The config file for InDel is in a BED format tabular separated plain text file, with coordinates started from 0. The REF and ALT are the bases for InDel in the reference genome and variant. An example file is given for 8-kb fragments containing human CODIS STR. For Details, please read our paper. There is one fixed headline and 5 columns in each line. Then each subsequent line is for one InDel site. If multiple InDel for the same coordinate, list one of them is enough. 

    Chrom	ChromStart	ChromEnd	REF	ALT
    chr1	230764747	230764748	GA	G
    chr1	230764749	230764750	T	TC
    chr2	1485704	1485705	C	CT
    chr2	1486795	1486803	TAGTGGGG	T
    chr2	1487522	1487526	CAGG	C
    chr2	1489924	1489975	GCACACAGGAGGAGTCACGACAGAGCAGTGTAAGAGCCGCCACGTGGGTCC	G

## Configure format for SNP
The config file for SNP is in a tabular separated plain text file, with coordinates started from 1.  There is one fixed headline and 5 columns in each line. Then each subsequent line is for one SNP site. An example file is given for 8-kb fragments containing 20 human CODIS STRs. the "ID" column is the name of a DNA fragment, which uses the same name as one of the STR in this 8-kb fragment as STR configure. All SNPs with the same "ID" will be put together with STR alleles in the same fragment once the ID of STR is shared between SNP and STR configure file. 

    #CHROM	 POS 	ID	REF	ALT
    chr1	230765010	D1S1656	A	G
    chr1	230765280	D1S1656	G	A
    chr1	230765336	D1S1656	C	T
    chr1	230765352	D1S1656	G	C
 
 ## Output of MacroHapCaller

Macrohaplotypes and supported read count in tabular text file .tsv. 
E.g., the two macrohaplotypes around FBI's CODIS loci D2S441 in 8 kb PacBio HiFi reads of benchmark reference hg002. Locus details are at [FBI](https://www.fbi.gov/how-we-can-help-you/dna-fingerprint-act-of-2005-expungement-policy/codis-and-ndis-fact-sheet)
     
The first line shows the total reads for this locus D3S1358. The position lines present the coordinate of each targeted variant site on chromosomes in the human genome Hg38. Then the macrohaplotypes, including three parts of SNPs, InDels and STRs separated by ";".

     #Total hapVar: 	D3S1358	47621	
     #Markername	Counts	HapVarLen	hapVar(s)
      #			  Position:45534683,45534822,45535067,45535094,45535770,45535776,45535998,45536129,45536613,45536683,45536800,45536841,45537042,45537158,45537183,45537457,45537777,45537899,45538041,45538207,45538590,45538748,45539142,45539154,45539315,45539413,45539452,45539580,45539887,45539923,45540056,45540467,45540750,45540754,45540758,45541155,45541805,45541806,45537166,45537617,45539736,45540709
      D3S1358	13843	108	T,A,T,C,A,C,T,A,G,A,T,A,A,G,G,C,C,C,T,T,C,G,T,G,A,C,G,G,T,A,G,C,A,A,A,T,C,G;TG,CA,TA;TCTATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA
      D3S1358	13600	101	G,G,T,T,A,A,T,G,C,G,C,A,A,T,G,T,G,T,T,T,C,A,T,G,A,T,A,A,T,A,G,C,G,G,A,T,C,G;T,C,T;TCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA

    
   ## Visualization 
   The macrohaplotypes can be visualized or edited in a text processor, Microsoft excels, and spreadsheet.
   
   The graphically view and compare the macrohaplotypes, please use our graphic tool: [USAT](https://github.com/XuewenWangUGA/USAT)
   
  ## Other tools for processing bam files
  These tools are designed to process bam in any computing platform, bam processing for windows. It can be a Windows version of some functions of `samtools`. Please make sure you have the latest version of JAVA standard environment installed (Java 17 or higher version).  This tool has been tested on Windows, Linux, and MacOS. 
  Java environment download link: https://www.oracle.com/java/technologies/downloads/. 

   #### Subset of a bam file:
  
for help:   `java -jar BamSubset.jar`
 
Usage: 
  
`java -jar BamSubset.jar    Bed_file     Input.bam     reduced_input.bam`
  
 
  #### Index and sort a bam file:

for help:

`java -jar BamAlignSortIndex.jar`  

Usage:

`java -jar BamAlignSortIndex.jar input.bam`


  
  ## Fund
  
  This work was supported in part by award 15PNIJ-21-GG-04159-RESS, awarded by the National Institute of Justice, Office of Justice Programs, U.S. Department of Justice.
  
   ## Citation
   coming soon






