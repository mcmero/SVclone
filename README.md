# README #

This package allows the clustering of subclonal structural variations. The package is divided into three components: process, filter and cluster. The process submodule counts the relevant reads from breakpoint locations, given a list of structural variants and a BAM file. The filter submodule filters the obtained read counts and prepares them for clustering; optionally, it can also prepare SNVs for clustering. The clustering submodule performs the clustering of the variants, optionally with SNVs.  

### How do I get set up? ###

Ensure you have the following dependencies installed:

Pre-process/Process modules

* [Numpy](http://www.numpy.org/) - install for python 2
* [PySam](http://pysam.readthedocs.org/en/latest/)
* [PyVCF](https://pyvcf.readthedocs.org/en/latest/)

Filter module

* [Numpy](http://www.numpy.org/) - install for python 2
* [Pandas](http://pandas.pydata.org/)
* [PyVCF](https://github.com/jamescasbon/PyVCF)

Cluster module

* [Numpy](http://www.numpy.org/) - install for python 2
* [Pandas](http://pandas.pydata.org/)
* [PyMC](https://github.com/pymc-devs/pymc)
* [ipython](https://pypi.python.org/pypi/ipython)
* [matplotlib](http://matplotlib.org/)

Install like so:

    python setup.py install

### Identify step (SV annotation) ###

If your SVs are in VCF or Socrates format, or are lacking direction or classification information, you will have to run them through the _preprocess_ step.  

    ./SVClone.py identify -i <svs> -b <indexed_bamfile> -o <output_base_name>

Input is expected in VCF format. Each defined SV must have a matching mate, given in the MATEID value in the INFO section. Input may also be entered in Socrates or simple format. Simple format is as follows:

```
bp1_chr	bp1_pos	bp2_chr	bp2_pos
22	18240676	22	18232335
22	19940482	22	19937820
22	21383572	22	21382745
22	21395573	22	21395746
```

Input MUST BE SORTED for Socrates and simple input methods. (bp1 < bp2 position and bp1s should be in chromosome/position order.)

Optionally a classification field may be specified with --sv_class, or directions for each break can be specified, if included in the input file, by specifying --use_dir. Note that the --use_dir flag does not work for VCF files (yet). 

#### Required Parameters ####

* -i or --input : structural variants input file (see above).
* -b or --bam : bam file with corresponding index file.
* -s or --sample : Sample name. Will create processed output file as <outdir>/<sample>_svinfo.txt, parameters output as <outdir>/<sample>_params.txt.

#### Optional Parameters ####

* -o or --outdir <outdir> : output directory to create files. Default: the sample name.
* -cgf or --config <config.ini>: SVClone configuration file with additional parameters (svclone_config.ini is the default).
* -r or --read_len : Read length of the bam file. Will be inferred if not specified (WARNING: if your read lengths are not constant, you will have to manually specify this parameter).
* --simple : Run using simple file format type as input (see File Formats section).
* --socrates : Use a Socrates-format style file as SV calls input (the input file must contain headers, these can be specified in the SVClone/SVProcess/parameters.py file).
* --use_dir : Whether to use breakpoint direction in the input file (where it must be supplied).
* --filter_repeats : Repeat types to filter out (can be a comma-separated list). SOCRATES INPUT ONLY.
* --sv_class_field : If your SV list has classifications and you would like to use them, specify the field name. 
* --min_mapq : Filter out SVs with lower average MAPQ than this value. SOCRATES INPUT ONLY (default 0).
* --trust_sc_pos : Use specified breaks without checking for differing soft-clip consensus position. Cannot be skipped if directionality must be inferred. If your SV caller offsets breaks due to micro-homology, e.g. Socrates, using this option is not recommended.
* --blacklist <file.bed> : Takes a list of intervals in BED format. Skip processing of any break-pairs where either SV break-end overlaps an interval specified in the supplied bed file. 

### Count step (SV read counting) ###

Run SV processing submodule to obtain read counts around breakpoints on each sample BAM file like so:

    ./SVClone.py process -i <svs> -b <indexed_bamfile> -o <output_base_name>

The classification strings are not used by the program, except for DNA-gain events (such as duplications). The classification names for these types of SVs should be specified in the svclone_config.ini file (see configuration file section).

#### Required Parameters ####

* -i or --input : structural variants input file. This should be the output file from the Identify step. 
* -b or --bam : bam file with corresponding index file.
* -s or --sample : Sample name. Will create processed output file as <outdir>/<sample>_svinfo.txt, parameters output as <outdir>/<sample>_params.txt.

#### Optional Parameters ####

* -o or --outdir <outdir> : output directory to create files. Default: the sample name.
* -cgf or --config <config.ini>: SVClone configuration file with additional parameters (svclone_config.ini is the default).
* -d or --depth <value> (default = 50) : average depth of coverage for corresponding BAM file. Used to skip areas of higher than expected coverage. Note that highly accurate figures are not necessary, a rough estimate will do.
* -r or --read_len <value> (automatically inferred if not supplied) :  specify if the read length is known, otherwise the program will infer this through the supplied bam file. If you have varying read sizes, we suggest you trim these reads to a constant size. The program will likely crash if it detects different read sizes, unless this parameter is supplied. 
* -v or --insert_mean <value> (automatically inferred if not supplied) : the average fragment or template length. Specify if known, otherwise the program will infer this through the first 50,000 reads in the supplied bam file.
* --insert_std <value> (automatically inferred if not supplied) : specify if the insert standard deviation is known, otherwise the program will infer this through the supplied bam file.
* --write_anomalous : Anomalous reads are reads that cross the SV boundary, but are not counted as supporting SV reads or normal reads. This flag will write a bam file with the reads and will recount anomalous reads accurately. Useful for diagnostic purposes. 

### Filter step (Filter SVs and/or SNVs and attach CNV states) ###

To filter the data obtained from the SV counting program and/or filter SNV data, can be done like so:

    ./SVClone.py filter -i <sv_info.txt> -s <sample_name> -c <battenberg_subclones.txt> --snvs <snvs_file> -o <output_directory>

Running the flat clustering approach on a single sample (the only currently supported method), can be done like so:

#### Required Parameters ####

* -s or --sample <name> : sample name, currently only a single sample is supported. WARNING: if clustering using mutect SNVs, the sample name must match the sample name in the vcf file.
* -i or --input <svinfo.txt> : sv info file from SV pre-processing script output.

#### Optional Parameters ####

Note that read length and insert sizes used by the filter step are provided as outputs from the pre-processing script (<out>_params.txt), based on the first 1000 sampled reads in the bam file. 

* -o or --outdir <outdir> : output directory to create files. Default: the sample name.
* -cgf or --config <config.ini>: SVClone configuration file with additional parameters (svclone_config.ini is the default).
* --params <params.txt> : Parameters file from processing step containing read information. If not supplied, the default search path is <outdir>/<sample>_params.txt'
* -c or --cnvs <cnv file> : Battenberg subclones.txt file containing segmented copy-numbers for patient sample. If not supplied, will assume that all regions are copy-number neutral (Major = 1, Minor = 1).
* -p <file> : Tumour purity and ploidy in text file. Purity must be between 0 and 1. Ploidy can be a floating point number.
* -g or --germline <germline_svinfo.txt> : Germline SVs; will filter out any tumour SVs which have >1 supporting read in the germline. Excpects the same input format as the sv info file (you can run sv_process on the tumour SVs against the germline bam file to obtain this).
* --neutral : Keep only copy-number neutral SVs.
* --snvs <snv_file> : SNVs in VCF format to (optionally) compare the clustering with SVs.
* --snv_format <sanger,mutect,mutect_callstats> (default = sanger) : Specify VCF input format (only if clustering SNVs).
* --minsplit <value> (default = 1) : Require at least N split reads to keep SV.
* --minspan <value> (default = 1) : Require at least N spanning reads to keep SV.
* --sizefilter <value> (default (readlen * 2) + insert_mean) : Filter out SVs below this size; the default minimum size is the fragment size. 
* --filter_outliers <value> : Filter out SVs with depth values that are considers outliers, based on the copy-number adjusted distribution of depths.
* --valid_chroms : Filter out SVs not on standard chromosomes (i.e. mapping to contigs). The list of valid chromosomes is specified in the SVClone/parameters.py file. 
* --min_depth : Filter any variants with total depth below this value (default = 4). Applies to both SVs and SNVs.
* --blacklist <file.bed> : Takes a list of intervals in BED format. Skip processing of any break-pairs where either SV break-end overlaps an interval specified in the supplied bed file.
* --strict_cnv_filt : Removes variants with no matched CNV state, otherwise assumes the CNV state is ploidy/2 for major and minor (when round(ploidy) < 2, state becomes 1-0).

### Purity/ploidy file ###

Purity/ploidy file must have the following layout (tab-separated):

```
sample	purity	ploidy
<sample>	<float between 0 - 1>	<positive float or integer>
```
### Clustering SVs ###

Once we have the filtered SV and/or SNV counts, we can run the clustering:

    ./SVClone.py cluster -s <sample_name> -o <outdir>

#### Required Parameters ####

* -s or --samples : sample names, currently only a single sample is supported.
* -o or --outdir : output directory (must be the same as the output directory from the filter step). 

#### Optional Parameters ####

Note that read length and insert sizes are provided as outputs from the pre-processing script (<out>_params.txt), based on the first 1000 sampled reads in the bam file. 

* --params <params.txt> : Parameters file from processing step containing read information. If not supplied, the default search path is <outdir>/<sample>_params.txt'
* -n or --n_runs <value> (default = 1) : number of times to run complete rounds of MCMC sampling (does not set the number of MCMC iterations, but the number of times the clustering runs are performed). Each run will have distinct results. 
* -t or --n_iter <value> (default = 10000) : the number of MCMC iterations to perform.
* --burn <value> (default = 0) : MCMC burnin parameter
* --thin <value> (default = 1) : MCMC thinning parameter
* --merge : whether to perform cluster merging.
* --map : use maximum a-posteriori fitting (may significantly increase runtime).
* --cocluster : cluster SVs and SNVs together.
* --no_adjust : do not adjust read counts based on different classes of SV events.

#### File Formats ####

If you run the program with the --simple argument, you can also provide SVs in a text-delimited format as follows:

```
bp1_chr	bp1_pos	bp1_dir	bp2_chr	bp2_pos	bp2_dir	classification
22	18240676		-	22	18232335	-	INV
22	19940482		-	22	19937820    +	DEL
22	21383572		-	22	21382745	+	DUP
22	21395573		+	22	21395746	+	INV 
```

Otherwise specify --socrates if your output is from the Socrates SV called (columns must be named, field names can be tweaked in the configuration file). If using the VCF format, SV pairs must be linked with the MATEID INFO field. 

#### Configuration file ####

If customisation is required for constants, these can be modified in the svclone_config.ini, or a new config file can be specified and each step run with the --config or -cfg flag. 

```
[GlobalParameters]
threshold: 6
# ^ "wobble length" tolerance threshold which we allow breaks to be inexact
germline_threshold: 10
# ^ bp threshold that a germline and tumour SV must match to be considered the same event
norm_overlap: 10   
# ^ minimum basepairs a "normal" read must overlap break to be counted
sc_len: 10
# ^ minimum basepairs a supporting read must be softclipped over the break
max_cn: 10
# ^ maximum possible copy-number
mean_cov: 50
# ^ mean coverage of the bam
support_adjust_factor: 0.25
# ^ base scaling factor for supporting reads = 1 + (support_adjust_factor * purity)
sv_offset: 100000
# ^ SVs offset by this amount of bp when matching CNVs

[ValidationParameters]
chroms: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y

[SVClasses]
# Naming conventions used to label SV types
inversion_class: INV
deletion_class: DEL
dna_gain_class: DUP,INTDUP
dna_loss_class: DEL,INV,TRX
itrx_class: INTRX

[SocratesFields]
# Column names used by Socrates SV caller (not needed if using different SV caller)
bp1_pos: C1_anchor
bp1_dir: C1_anchor_dir
bp2_pos: C1_realign
bp2_dir: C1_realign_dir
avg_mapq1: C1_avg_realign_mapq
avg_mapq2: C2_avg_realign_mapq
repeat1: repeat1
repeat2: repeat2

[ClusterParameters]
phi_limit: 2
# ^ max CCF
clus_limit: 25
# ^ max number of possible clusters
subclone_diff: 0.10
# ^ max difference in CCF between subclones to not recluster
hpd_alpha: 0.05
# ^ credible interval (for computing highest posterior density interval)

[BetaParameters]
# Control cluster sensitivity
alpha: 0.9
beta: 1

```

### Post Processing ###

To create some helpful plots for SV clustering results (currently only SV output is supported), use the following Rscript:

```
Rscript post_process_sv_only.R <workingdir> <sample_id> <run_number>
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

For queries, contact cmerom[at]gmail.com
