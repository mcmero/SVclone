# README #

The script obtains all information necessary from a bam file and a list of structural variations in order to accurately count variant allele frequeny (VAF). 

### How do I get set up? ###

Ensure you have the following dependencies installed:

* [Numpy](http://www.numpy.org/) - install for python 2
* [PySam](http://pysam.readthedocs.org/en/latest/)

Install like so:

    python setup.py install

### Pre-processing of SVs ###

Run SV pre-processing on each sample BAM file like so:

    python -m sv_process.cmd -i <svs.txt> -b <indexed bamfile> -o <output base name> -d <average coverage>

#### Required Parameters ####

* -i or --input : structural variants input file. See file formats section for more details.
* -b or --bam : bam file with corresponding index file.
* -o or --out : output base name. Will create processed output file as <name>_svinfo.txt, parameters output as <name>_params.txt and database output as <name>_svinfo.db
* -d or --depth (floating-point value) : average depth of coverage for corresponding BAM file. Used to skip areas of higher than expected coverage.

#### Optional Parameters ####

* -sc or --softclip (default = 25) : reads must overlap by this many basepairs to be counted as supporting the break, or being a non-supporting normal read lying across the break.
* -cn or --max_cn (default = 15) : maximum expected copy-number. Will skip any areas with more average depth higher than <depth> * <max_cn>
* --read_len (inferred automatically by default) :  specify if the read length is known, otherwise the program will infer this through the supplied bam file.
* --insert_mean (inferred automatically by default) : specify if the insert length is known, otherwise the program will infer this through the supplied bam file.
* --insert_std (inferred automatically by default) : specify if the insert standard deviation is known, otherwise the program will infer this through the supplied bam file.
* 

#### Beyond Advanced Parameters ####

The package also contains a parameters.py file which has the following hard-coded parameters. Modify these with care.
```
tr      = 5    # if soft-clipped by less than these number of bases at this end, is not a "true" soft-clip
window  = 500  # base-pair window considered to the left and right of the break when processing reads
```
#### File Formats ####

Input is expected in VCF format as follows:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOUR
1	44496801	3_1	C	C]11:574949]	85	.	SVTYPE=BND;MATEID=3_2;IMPRECISE;CIPOS=0,2;CIEND=2,0;HOMSEQ=AC;HOMLE`N=2;TRDS=C26GFACXX130520:7:1214:2378:56695/1,D282HACXX130520:1:1212:8778:24496/1,C26GFACXX130520:7:1212:2158:48685/2,D282HACXX130520:2:1104:8301:15718/2,D282HACXX130520:1:2114:19802:53112/1;GENE=SLC6A9;SID=SLC6A9;TID=ENST00000372310;AS=-1;RGN=5UTRintron;RGNNO=1;RGNC=14;	RC	0	5
1	66558497	4_1	A	A[1:66974085[	100	.	SVTYPE=BND;MATEID=4_2;IMPRECISE;CIPOS=0,1;CIEND=0,1;HOMSEQ=A;HOMLEN=1;TRDS=D282HACXX130520:1:1108:16793:100498/2,D282HACXX130520:2:1314:4766:40399/1,C26GFACXX130520:7:1310:15450:41573/2;GENE=PDE4B;SID=PDE4B;TID=ENST00000329654;AS=1;EPH=2;RGN=intron;RGNNO=3;RGNC=17;	RC	0	3
1	66974085	4_2	A	]1:66558497]A	100	.	SVTYPE=BND;MATEID=4_1;IMPRECISE;CIPOS=0,1;CIEND=0,1;HOMSEQ=A;HOMLEN=1;TRDS=D282HACXX130520:1:1108:16793:100498/2,D282HACXX130520:2:1314:4766:40399/1,C26GFACXX130520:7:1310:15450:41573/2;	RC	0	3
1	150061074	6_1	C	[1:151855796[TTGAAC	100	.	SVTYPE=BND;MATEID=6_2;FFV=700;TRDS=D282HACXX130520:2:2111:7011:83262/2,C26GFACXX130520:7:2103:6951:62396/2,C26GFACXX130520:7:1313:15644:77031/1,D282HACXX130520:1:2311:7590:76009/2;GENE=VPS45;SID=VPS45;TID=ENST00000369130;AS=1;RGN=intron;RGNNO=10;RGNC=15;	RC	0	4
1	151855796	6_2	G	[1:150061074[TTGAAG	100	.	SVTYPE=BND;MATEID=6_1;FFV=700;TRDS=D282HACXX130520:2:2111:7011:83262/2,C26GFACXX130520:7:2103:6951:62396/2,C26GFACXX130520:7:1313:15644:77031/1,D282HACXX130520:1:2311:7590:76009/2;GENE=THEM4;SID=THEM4;TID=ENST00000368814;AS=-1;PH=2;RGN=intron;RGNNO=4;RGNC=6;	RC	0	4
```

Each defined SV must have a matching mate, given in the MATEID value in the INFO section.

Optionally, if you run the program with the --simple argument, you can also provide SVs in a simple text-delimited format as follows:

```
bp1_chr	bp1_pos	bp2_chr	bp2_pos	classification
22	18240676	22	18232335	INV
22	19940482	22	19937820    INV
22	21383572	22	21382745	INV
22	21383573	22	21382746	INV 
```

### Calculating Coverage ###

Coverage can be approximately calculated using a script such as:
```
#!/bin/sh
bam=$1
genome=$2

bedtools random -g $genome -n 1000 -l 1000 > rand_intervals_1kb.bed
bedtools sort -chrThenSizeA -i rand_intervals_1kb.bed > rand_intervals_1kb.sort.bed
coverageBed -abam $bam -b rand_intervals_1kb.sort.bed > cov.txt

#in coverage_script.R
rlen = 100
interval = 1000
x <- read.delim(‘cov.txt’,header=F,stringsAsFactors=F)
print((mean(x[x$V3!=0,’V3’])*rlen)/interval
```
### Who do I talk to? ###

For queries, contact cmerom@gmail.com
