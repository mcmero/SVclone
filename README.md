# README #

This package allows the clustering of structural variation cancer cell fractions (CCFs). The package is divided into four components: annotate, count, filter and cluster. The annotate step infers directionality of each breakpoint (if not supplied), recalibrates breakpoint position to the soft-clip boundary and subsequently classifies SVs using a rule-based approach. The count step counts the variant and non-variant reads from breakpoint locations. Both the annotate and count steps utilise BAM-level information. The filter step removes SVs based on a number of adjustable parameters and prepares the variants for clustering. SNVs can also be added at this step as well as CNV information, which is matched to SV and SNV loci. Post-processing scripts are also included to aid in visualising the clustering results. 

### How do I get set up? ###

Ensure you have the following dependencies installed:

Annotate/count steps:

* [Numpy](http://www.numpy.org/) - install for python 2
* [PySam](http://pysam.readthedocs.org/en/latest/)
* [PyVCF](https://pyvcf.readthedocs.org/en/latest/)

Cluster module

* [Numpy](http://www.numpy.org/) - install for python 2
* [Pandas](http://pandas.pydata.org/)
* [PyMC](https://github.com/pymc-devs/pymc)
* [ipython](https://pypi.python.org/pypi/ipython)
* [matplotlib](http://matplotlib.org/)
* [SciPy](https://http://www.scipy.org/)

Install like so:

    python setup.py install

Adding --user if a local install is required. 

### Example data ###

Example data is provided to test your SVclone installation. Run as:

    ./run_example.sh

### Annotate step ###

An indexed whole-genome sequencing BAM and a list of paired breakpoints from an SV caller of choice is required. This step is required for clustering of SVs, however, classifiation and directionality information from your choice of SV caller can be used rather than being inferred. 

    ./SVclone.py annotate -i <sv_input> -b <indexed_bamfile> -s <sample_name>

Input is expected in VCF format (directionality inferred from the ALT field is also supported). Each defined SV must have a matching mate, given in the MATEID value in the INFO section. Input may also be entered in Socrates or simple format (must be specified with --sv_format simple or socrates). Simple format is as follows:

```
chr1	pos1	chr2	pos2
22	18240676	22	18232335
22	19940482	22	19937820
22	21383572	22	21382745
22	21395573	22	21395746
```

or, if directionality or classification information is available:

```
chr1	pos1	dir1	chr2	pos2	dir2	classification
22	18240676		-	22	18232335	-	INV
22	19940482		-	22	19937820    +	DEL
22	21383572		-	22	21382745	+	DUP
22	21395573		+	22	21395746	+	INV
```

Input MUST BE SORTED for Socrates and simple input methods. (bp1 < bp2 position and bp1s should be in chromosome/position order.)

Optionally a classification field may be specified in the 'sv_class_field' parameter in the configuration file, additionally, to specify directionality set the parameter 'use_dir' to True. A blacklist can also be supplied at this step to not process areas to remove SVs where any of its breakpoints fall into one of these areas. 

#### Required Parameters ####

* -i or --input : structural variants input file (see above).
* -b or --bam : bam file with corresponding index file.
* -s or --sample : Sample name. Will create processed output file as <out>/<sample>_svinfo.txt, parameters output as <out>/<sample>_params.txt.

#### Optional Parameters ####

* -o or --out <directory> : output directory to create files. Default: the sample name.
* -cgf or --config <config.ini> : SVClone configuration file with additional parameters (svclone_config.ini is the default).
* --sv_format <vcf, simple, socrates> : input format of SV calls, VCF by default, but may also be simple (see above) or from the SV caller Socrates. 
* --blacklist <file.bed> : Takes a list of intervals in BED format. Skip processing of any break-pairs where either SV break-end overlaps an interval specified in the supplied bed file. Using something like the [DAC blacklist](https://www.encodeproject.org/annotations/ENCSR636HFF/) is recommended.

### Count step ###

Run SV processing submodule to obtain read counts around breakpoints on each sample BAM file like so:

    ./SVClone.py count -i <svs> -b <indexed_bamfile> -s <sample_name>

The classification strings are not used by the program, except for DNA-gain events (such as duplications). The classification names for these types of SVs should be specified in the svclone_config.ini file (see configuration file section).

#### Required Parameters ####

* -i or --input : structural variants input file. This should be the output file from the annotate step.
* -b or --bam : bam file with corresponding index file.
* -s or --sample : Sample name. Will create processed output file as <out>/<sample>_svinfo.txt, parameters output as <out>/<sample>_params.txt.

#### Optional Parameters ####

* -o or --out <directory> : output directory to create files. Default: the sample name.
* -cgf or --config <config.ini>: SVClone configuration file with additional parameters (svclone_config.ini is the default).

### Filter step (Filter SVs and/or SNVs and attach CNV states) ###

To filter the data obtained from the SV counting program and/or filter SNV data, can be done like so:

    ./SVClone.py filter -i <sv_info.txt> -s <sample_name>

Note that read length and insert sizes used by the filter step are provided as outputs from the count step (<out>/read_params.txt), based on the first 50,000 sampled reads in the bam file.

#### Required Parameters ####

* -s or --sample <name> : sample name, currently only a single sample is supported. WARNING: if clustering using mutect SNVs, the sample name must match the sample name in the vcf file.
* -i or --input <svinfo.txt> : sv info file from SV count step.

#### Optional Parameters ####

* -o or --out <out> : output directory to create files. Default: the sample name.
* -cgf or --config <config.ini>: SVClone configuration file with additional parameters (svclone_config.ini is the default).
* --params <params.txt> : Parameters file from processing step containing read information. If not supplied, the default search path is <out>/<sample>_params.txt'
* -c or --cnvs <cnv file> : Battenberg subclones.txt file containing segmented copy-numbers for patient sample. If not supplied, will assume that all regions are copy-number neutral (Major = 1, Minor = 1).
* -p <file> or --purity_ploidy <file>: Tumour purity and ploidy in tab-separated text file. Purity must be between 0 and 1. Ploidy can be a floating point number. Column names must be labelled 'sample', 'purity' and 'ploidy' (without quotes). Row 2 must contain the sample name and purity and ploidy values respectively.
* -g or --germline <germline_svinfo.txt> : Germline SVs; will filter out any tumour SVs which have >1 supporting read in the germline. Excpects the same input format as the sv info file (you can run sv_process on the tumour SVs against the germline bam file to obtain this).
* --snvs <snv_file> : SNVs in VCF format to (optionally) compare the clustering with SVs.
* --snv_format <sanger, mutect, mutect_callstats> (default = sanger) : Specify VCF input format (only if clustering SNVs).
* --blacklist <file.bed> : Takes a list of intervals in BED format. Skip processing of any break-pairs where either SV break-end overlaps an interval specified in the supplied bed file. Using something like the [DAC blacklist](https://www.encodeproject.org/annotations/ENCSR636HFF/) is recommended.
* --subsample <integer> : (SNVs only); subsample N variants from total filtered input.
* --seed <integer> : (only valid if using subsampling) integer seed, can be set to ensure subsampling is replicable.

### Purity/ploidy file ###

Purity/ploidy file must have the following layout (tab-separated):

```
sample	purity	ploidy
<sample>	<float between 0 - 1>	<positive float or integer>
```
### Clustering SVs ###

Once we have the filtered SV and/or SNV counts, we can run the clustering:

    ./SVClone.py cluster -s <sample_name>

#### Required Parameters ####

* -s or --samples : sample names, currently only a single sample is supported.

#### Optional Parameters ####

* -o or --out : output directory (sample name by default)
* -cgf or --config <config.ini>: SVClone configuration file with additional parameters (svclone_config.ini is the default).
* --params <params.txt> : Parameters file from processing step containing read information. If not supplied, the default search path is <out>/<sample>_params.txt'
* --snvs <filtered_snvs> : SNVs output from filter step (automatically detected if present). 
* --map : use maximum a-posteriori fitting (may significantly increase runtime).
* --cocluster : cluster SVs and SNVs together.
* --no_adjust : do not adjust read counts based on different classes of SV events.
* --seeds <one integer per run, comma separated> : set a random seed per run in order to replicate pyMC runs. 
* --XX and --XY : overwrite the config file genotype with XX or XY. 

#### Configuration file ####

If customisation is required for parameters, these can be modified in the svclone_config.ini, or a new config file can be specified and each step run with the --config or -cfg flag for all steps.
