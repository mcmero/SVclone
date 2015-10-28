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

### Pre-processing of SVs ###

If your SVs are in VCF or Socrates format, or are lacking direction or classification information, you will have to run them through the _preprocess_ step.  

    ./SVClone.py preprocess -i <svs> -b <indexed_bamfile> -o <output_base_name>

Input is expected in VCF format. Each defined SV must have a matching mate, given in the MATEID value in the INFO section. Input may also be entered in Socrates or simple format. Simple format is as follows:

```
bp1_chr	bp1_pos	bp2_chr	bp2_pos
22	18240676	22	18232335
22	19940482	22	19937820
22	21383572	22	21382745
22	21395573	22	21395746
```

Input MUST BE SORTED for Socrates and simple input methods. (bp1 < bp2 position and bp1s should be in chromosome/position order.)

Optionally a classification field may be specified with --sv_class, or directions for each break can be specified, if included in the input file, by specifying --use_dir. 

#### Required Parameters ####

* -i or --input : structural variants input file (see above).
* -b or --bam : bam file with corresponding index file.
* -o or --out : output base name. Will create processed output file as <name>_svinfo.txt, parameters output as <name>_params.txt and database output as <name>_svinfo.db

#### Optional Parameters ####

* -md or --max_dep <value> (default = 1000) : maximum depth to consider when extracting direction information. Will skip all locations where there are more than this many reads at the locus. 
* --simple : run using simple file format type as input (see File Formats section).
* --socrates : use a Socrates-format style file as SV calls input (the input file must contain headers, these can be specified in the SVClone/SVProcess/parameters.py file).
* --use_dir : whether to use breakpoint direction in the input file (where it must be supplied).
* --filter_repeats : Repeat types to filter out (can be a comma-separated list). SOCRATES INPUT ONLY.
* --sv_class_field : If your SV list has classifications and you would like to use them, specify the field name. 
* --min_mapq : Filter out SVs with lower average MAPQ than this value. SOCRATES INPUT ONLY (default 0).


### Processing of SVs ###

Run SV processing submodule to obtain read counts around breakpoints on each sample BAM file like so:

    ./SVClone.py process -i <svs> -b <indexed_bamfile> -o <output_base_name>

The process step expects input in the following format (output from preprocess):

```
bp1_chr	bp1_pos	bp1_dir	bp2_chr	bp2_pos	bp2_dir	classification
22	18240676		-	22	18232335	-	INV
22	19940482		-	22	19937820    +	DEL
22	21383572		-	22	21382745	+	DUP
22	21395573		+	22	21395746	+	INV 
```

The classification strings are not used by the program, except for DNA-gain events (such as duplications). The classification names for these types of SVs should be specified in the SVClone/parameters.py file: 

```
dna_gain_class = ['DUP','INTDUP']
```

#### Required Parameters ####

* -i or --input : structural variants input file. See file formats section for more details.
* -b or --bam : bam file with corresponding index file.
* -o or --out : output base name. Will create processed output file as <name>_svinfo.txt, parameters output as <name>_params.txt and database output as <name>_svinfo.db

#### Optional Parameters ####

* -d or --depth <value> (default = 50) : average depth of coverage for corresponding BAM file. Used to skip areas of higher than expected coverage. Note that highly accurate figures are not necessary, a rough estimate will do.
* -sc or --softclip <value> (default = 25) : reads must overlap by this many basepairs to be counted as supporting the break, or being a non-supporting normal read lying across the break.
* -cn or --max_cn <value> (default = 15) : maximum expected copy-number. Will skip any areas with more average depth higher than <depth> * <max_cn>
* -r or --read_len <value> (automatically inferred if not supplied) :  specify if the read length is known, otherwise the program will infer this through the supplied bam file. If you have varying read sizes, we suggest you trim these reads to a constant size. The program will likely crash if it detects different read sizes, unless this parameter is supplied. 
* -v or --insert_mean <value> (automatically inferred if not supplied) : the average fragment or template length. Specify if known, otherwise the program will infer this through the first 50,000 reads in the supplied bam file.
* --insert_std <value> (automatically inferred if not supplied) : specify if the insert standard deviation is known, otherwise the program will infer this through the supplied bam file.
* --simple : run using simple file format type as input (see File Formats section).
* --socrates : use a Socrates-format style file as SV calls input (the input file must contain headers).
* --use_dir : whether to use breakpoint direction in the input file (where it must be supplied).
* --filter_repeats : Repeat types to filter out (can be a comma-separated list). SOCRATES INPUT ONLY.
* --sv_class_field : If your SV list has classifications and you would like to use them, specify the field name. 

#### Beyond Advanced Parameters ####

The package also contains a parameters.py file which has the following hard-coded parameters. Modify these with care.

```
tr      = 6    # if soft-clipped by less than these number of bases at this end, is not a "true" soft-clip
window  = 500  # base-pair window considered to the left and right of the break when processing reads
```

#### File Formats ####


Optionally, if you run the program with the --simple argument, you can also provide SVs in a simple text-delimited format as follows:


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

### Filtering SVs ###

To filter the data obtained from the SV counting program and/or filter SNV data, can be done like so:

    ./SVClone.py filter -i <sv_info.txt> -s <sample_name> --params <params.txt> -c <battenberg_subclones.txt> --snvs <snvs_file> -o <output_directory>

Running the flat clustering approach on a single sample (the only currently supported method), can be done like so:

#### Required Parameters ####

* -s or --samples <name> : sample name, currently only a single sample is supported.
* -i or --input <svinfo.txt> : sv info file from SV pre-processing script output.
* -o or --outdir <outdir> : output directory to create files.

#### Optional Parameters ####

Note that read length and insert sizes used by the filter step are provided as outputs from the pre-processing script (<out>_params.txt), based on the first 1000 sampled reads in the bam file. 

* --params <params.txt> : Parameters file from processing step containing read information. If not supplied, the default search path is <outdir>/<sample>_params.txt'
* -c or --cnvs <cnv file> : Battenberg subclones.txt file containing segmented copy-numbers for patient sample. If not supplied, will assume that all regions are copy-number neutral (Major = 1, Minor = 1).
* -p or --purity <value> (default = 1.0) : tumour purity. A floating point number between 0 - 1 indicating the percentage of tumour cells.
* -g or --germline <germline_svinfo.txt> : Germline SVs; will filter out any tumour SVs which have >1 supporting read in the germline. Excpects the same input format as the sv info file (you can run sv_process on the tumour SVs against the germline bam file to obtain this).
* -y or --ploidy <value> (default = 1.0) : tumour ploidy, affects the phi parameter (base model assumes variants occur on one copy only). Setting the ploidy to 1 gives the raw frequencies (it ignores ploidy).
* --neutral : Keep only copy-number neutral SVs.
* --snvs <snv_file> : SNVs in VCF format to (optionally) compare the clustering with SVs.
* --snv_format <sanger,mutect,mutect_callstats> (default = sanger) : Specify VCF input format (only if clustering SNVs).
* --minsplit <value> (default = 1) : Require at least N split reads to keep SV.
* --minspan <value> (default = 1) : Require at least N spanning reads to keep SV.
* --sizefilter <value> (default (readlen * 2) + insert_mean) : Filter out SVs below this size; the default minimum size is the fragment size. 
* --filter_outliers <value> : Filter out SVs with depth values that are considers outliers, based on the copy-number adjusted distribution of depths.
* --valid_chroms : Filter out SVs not on standard chromosomes (i.e. mapping to contigs). The list of valid chromosomes is specified in the SVClone/parameters.py file. 

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
* --burn (default = 0) : MCMC burnin parameter
* --thin (default = 1) : MCMC thinning parameter
* --beta (default = "0.01,1") : lower and upper bounds respectively of uniform distribution of beta parameter used in Dirichlet Process model. 
* --merge : whether to perform cluster merging.
* --map : use maximum a-posteriori fitting (may significantly increase runtime).
* --cocluster : cluster SVs and SNVs together.
* --no_adjust : do not adjust read counts based on different classes of SV events.

#### Beyond Advanced Usage ####

The package also contains a parameters.py file which has the following hard-coded parameters. Modify these with care.

```
subclone_threshold      = 0.05 # throw out any subclones with frequency lower than this value
subclone_sv_prop        = 0.08 # remove any cluster groups with fewer than this proportion of SVs clustering together
subclone_diff           = 0.10 # merge any clusters within this range
clus_limit              = 20 # maximum number of clusters generated by dirichlet process
```

### Who do I talk to? ###

For queries, contact cmerom@gmail.com
