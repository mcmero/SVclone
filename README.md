# README #

This package is used to cluster structural variants of similar cancer cell fraction (CCF). SVclone is divided into five components: annotate, count, filter, cluster and post-assign. The annotate step infers directionality of each breakpoint (if not supplied), recalibrates breakpoint position to the soft-clip boundary and subsequently classifies SVs using a rule-based approach. The count step counts the variant and non-variant reads from breakpoint locations. Both the annotate and count steps utilise BAM-level information. The filter step removes SVs based on a number of adjustable parameters and prepares the variants for clustering. SNVs can also be added at this step as well as CNV information, which is matched to SV and SNV loci. Any variants that were filtered out, or left out due to sub-sampling can be added back using the post-assign step, which assigns each variant (which contains a >0 VAF and matching copy-number state, at minimum) to the most likely cluster (obtained from the cluster step). Post-processing scripts are also included to aid in visualising the clustering results.

### How do I get set up? ###

Assuming [Anaconda2](https://www.continuum.io/downloads) (or [Python 2.7.\*](https://www.python.org/downloads/) with [Numpy](http://www.numpy.org/), [SciPy](https://http://www.scipy.org/) and [PyMC](https://pymc-devs.github.io/pymc/INSTALL.html)) are installed, a quick install can be run as follows:

    curl -O https://github.com/mcmero/SVclone/archive/0.2.1.tar.gz
    tar -xvzf 0.2.1.tar.gz

    cd SVclone
    pip install -r requirements.txt
    python setup.py install

Check the INSTALL.md file for full instructions.

### Example data ###

Example data is provided to test your SVclone installation (data contains simulated clonal deletions). Run as:

    ./run_example.sh

NOTE: installing [R](https://www.r-project.org/) is required to run the post-processing script. See INSTALL.md for package dependencies.

You can check the following output plots:

    * tumour_p80\_DEL/tumour_p80_DEL_run_summary.pdf
    * tumour_p80_DEL/tumour_p80_DEL_run\*.pdf
    * tumour_p80_DEL/tumour_p80_DEL_best_run_svs_post_assign_best_fit.pdf

You can also test the simulated SNV data by running (this will take longer than running SVs only): 

    ./run_coclus_example.sh

The simulated data contains a 100% CCF clone and a 30% subclone.

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

Optionally a classification field may be specified in the 'sv_class_field' parameter in the configuration file, additionally, to specify directionality set the parameter 'use_dir' to True. For example:
```
chr1	pos1	dir1	chr2	pos2	dir2	classification
22	18240676		-	22	18232335	-	INV
22	19940482		-	22	19937820    +	DEL
22	21383572		-	22	21382745	+	DUP
22	21395573		+	22	21395746	+	INV
```

A blacklist can also be supplied at this step to not process areas to remove SVs where any of its breakpoints fall into one of these areas.

#### Annotate Output ####

The above input example also corresponds with the output of this step (output to <out>/<sample>_svin.txt), with an added SV ID present in the column. Events that are considered part of the same event will have the same ID (which may be multiple breakpoints).

#### Required Parameters ####

* -i or --input : structural variants input file (see above).
* -b or --bam : bam file with corresponding index file.
* -s or --sample : Sample name. Will create processed output file as <out>/<sample>_svinfo.txt, parameters output as <out>/<sample>_params.txt.

#### Optional Parameters ####

* -o or --out <directory> : output directory to create files. Default: the sample name.
* -cgf or --config <config.ini> : SVclone configuration file with additional parameters (svclone_config.ini is the default).
* --sv_format <vcf, simple, socrates> : input format of SV calls, VCF by default, but may also be simple (see above) or from the SV caller Socrates.
* --blacklist <file.bed> : Takes a list of intervals in BED format. Skip processing of any break-pairs where either SV break-end overlaps an interval specified in the supplied bed file. Using something like the [DAC blacklist](https://www.encodeproject.org/annotations/ENCSR636HFF/) is recommended.

### Count step ###

Run SV processing submodule to obtain read counts around breakpoints on each sample BAM file like so:

    ./SVclone.py count -i <svs> -b <indexed_bamfile> -s <sample_name>

The classification strings are not used by the program, except for DNA-gain events (such as duplications). The classification names for these types of SVs should be specified in the svclone_config.ini file (see configuration file section).

#### Count output ####

The count step will create a tab-separated <out>/<sample>_svinfo.txt file containing count information. For example:

```
ID	chr1	pos1	dir1	chr2	pos2	dir2	classification	split_norm1	norm_olap_bp1	span_norm1	win_norm1	split1	sc_bases1	total_reads1	split_norm2	norm_olap_bp2	span_norm2	win_norm2	split2	sc_bases2	total_reads2	anomalous	spanning	norm1	norm2	support	vaf1	vaf2
1	12	227543	+	12	228250	-	DEL	12	405	13	96	4	215	189	15	473	8	94	4	149	190	32	4	25	23	12	0.32432432432432434	0.34285714285714286
2	12	333589	+	12	338298	-	DEL	19	585	23	132	1	69	222	19	492	13	100	8	385	213	18	12	42	32	21	0.33333333333333331	0.39622641509433965
3	12	461142	+	12	465988	-	DEL	14	490	12	120	6	247	202	12	374	16	104	6	149	214	20	6	26	28	18	0.40909090909090912	0.39130434782608697
4	12	526623	+	12	554937	-	DEL	11	322	18	112	8	312	220	17	567	15	106	8	232	205	12	9	29	32	25	0.46296296296296297	0.43859649122807015
5	12	693710	+	12	696907	-	DEL	13	433	15	104	9	329	212	16	446	21	138	5	245	229	20	9	28	37	23	0.45098039215686275	0.38333333333333336
```

The output fields are briefly described:
* split: split read count at each locus
* split_norm/span_norm: number of normal split and spanning reads crossing the boundary at locus 1 and 2 respectively.
* norm_olap_bp: count of normal read base-pairs overlapping the break (for normal reads that cross the break boundary).
* win_norm: normal read count (no soft-clips, normal insert size) for all normal reads extracted from the locus window (+/- insert size from locus).
* sc_bases: count of soft-clipped bases corresponding to split reads crossing the break.
* norm: normal read count at each locus.
* spanning: number of spanning reads supporting the break.
* support: split1 + split2 + spanning.
* anomalous: reads not counted in any other category.
* vaf: support / (norm + support).

#### Required Parameters ####

* -i or --input : structural variants input file. This should be the output file from the annotate step.
* -b or --bam : bam file with corresponding index file.
* -s or --sample : Sample name. Will create processed output file as <out>/<sample>_svinfo.txt, parameters output as <out>/<sample>_params.txt.

#### Optional Parameters ####

* -o or --out <directory> : output directory to create files. Default: the sample name.
* -cgf or --config <config.ini>: SVclone configuration file with additional parameters (svclone_config.ini is the default).

### Filter step (Filter SVs and/or SNVs and attach CNV states) ###

To filter the data obtained from the SV counting program and/or filter SNV data, can be done like so:

    ./SVclone.py filter -i <sv_info.txt> -s <sample_name>

Note that read length and insert sizes used by the filter step are provided as outputs from the count step (<out>/read_params.txt), based on the first 50,000 sampled reads in the bam file.

#### Filter output ####

The filter step outputs the file <out>/<sample>_filtered_svs.tsv and/or <out>/<sample>_filtered_snvs.tsv depending on input. For SVs, the output is akin to the _svinfo.txt file format with added fields:

* norm_mean: average of norm1 and norm2
* gtype: copy-number state at the locus: "major, minor, CNV fraction" for example, "1,1,1.0". May be subclonal if battenberg input is supplied e.g. "1,1,0.7|2,1,0.3".
* adjusted_norm: selected normal locus with adjusted normal read counts (in case of DNA-gain).
* adjusted_support: total adjusted supporting read count.
* adjusted_depth: adjusted_norm + adjusted_support
* preferred_side: which side the support/CNV state comes from
* raw_mean_vaf: support / (mean(norm1, norm2) + support)
* adjusted_vaf: adjusted_support / adjusted_depth

For SNVs, example output looks like:

```
chrom	pos	gtype	ref	var
1	44395	1,1,1.0	33.0	15.0
1	4865339	1,1,1.0	23.0	25.0
1	13846233	1,1,1.0	28.0	25.0
1	33976797	1,1,1.0	30.0	19.0
1	51346133	1,1,1.0	33.0	22.0
```

Where:
* ref: total reference allele reads at locus.
* var: total variant allele reads at locus.

#### Required Parameters ####

* -s or --sample <name> : sample name, currently only a single sample is supported. WARNING: if clustering using mutect SNVs, the sample name must match the sample name in the vcf file.
* -i or --input <svinfo.txt> : sv info file from SV count step.

#### Optional Parameters ####

* -o or --out <out> : output directory to create files. Default: the sample name.
* -cgf or --config <config.ini>: SVclone configuration file with additional parameters (svclone_config.ini is the default).
* --params <params.txt> : Parameters file from processing step containing read information. If not supplied, the default search path is <out>/<sample>_params.txt'
* -c or --cnvs <cnv file> : Battenberg subclones.txt file containing segmented copy-numbers for patient sample. If not supplied, will assume that all regions are copy-number neutral (Major = 1, Minor = 1).
* -p <file> or --purity_ploidy <file>: Tumour purity and ploidy in tab-separated text file. Purity must be between 0 and 1. Ploidy can be a floating point number. Column names must be labelled 'sample', 'purity' and 'ploidy' (without quotes). Row 2 must contain the sample name and purity and ploidy values respectively.
* -g or --germline <germline_svinfo.txt> : Germline SVs; will filter out any tumour SVs which have >1 supporting read in the germline. Excpects the same input format as the sv info file (you can run sv_process on the tumour SVs against the germline bam file to obtain this).
* --snvs <snv_file> : SNVs in VCF format to (optionally) compare the clustering with SVs.
* --snv_format <sanger, mutect, mutect_callstats> (default = sanger) : Specify VCF input format (only if clustering SNVs).
* --blacklist <file.bed> : Takes a list of intervals in BED format. Skip processing of any break-pairs where either SV break-end overlaps an interval specified in the supplied bed file. Using something like the [DAC blacklist](https://www.encodeproject.org/annotations/ENCSR636HFF/) is recommended.

### Purity/ploidy file ###

Purity/ploidy file must have the following layout (tab-separated):

```
sample	purity	ploidy
<sample>	<float between 0 - 1>	<positive float or integer>
```
### Clustering SVs ###

Once we have the filtered SV and/or SNV counts, we can run the clustering:

    ./SVclone.py cluster -s <sample_name>

#### Cluster output ####

SVclone creates output based on the PCAWG output specification. This includes (per run):

* number_of_clusters: number of clusters found.
* <sample>_copynumber.txt: each variant's copy-number state, including total copy-number and number of chromosomes (alleles) bearing the mutation.
* <sample>_multiplicity.txt: the total copy-number, the number of copies the variant occurs on, the different multiplicity options and the probability of each.
* <sample>_assignment_probability_table.txt: probability of each variant's assignment to each cluster, based on number of times the proportion that a variant occurs in a particular cluster over all MCMC iterations.
* <sample>_cluster_certainty.txt: each variant's most likely assignment to a particular cluster and corresponding average proportion (CCF x purity).
* <sample>_fit.txt: each IC metric's score, plus some extra metrics from PyMC (lnL, logp etc.).
* <sample>_subclonal_structure: clusters found, the number of variants per cluster, the proportion and CCF.

And a few more files unique to SVclone:

* <sample>_vaf_ccf.txt: variants with raw mean VAF, adjusted VAF (see filter step output), variant CCF derived from the trace, the transformed variant CCF and the cluster mean CCF/proportion.
* <sample>_most_likely_copynumbers.txt: contains the output from the PCAWG copynumber output format, plus the variants gtypes, pv (probability of sampling a variant read) and the pv deviance from VAF (high deviance suggests low CCF confidence for the variant).
* phi_trace.txt.gz: dump of the phi trace
* z_trace.txt.gz: dump of the z trace.
* alpha_trace.txt.gz: dump of the alpha trace (if alpha is not fixed).
* cluster_trace.png: a plot of the z, phi and alpha traces and a VAF histogram.

#### Required Parameters ####

* -s or --samples : sample names, currently only a single sample is supported.

#### Optional Parameters ####

* -o or --out : output directory (sample name by default)
* -cgf or --config <config.ini>: SVclone configuration file with additional parameters (svclone_config.ini is the default).
* --params <params.txt> : Parameters file from processing step containing read information. If not supplied, the default search path is <out>/<sample>_params.txt'
* --snvs <filtered_snvs> : SNVs output from filter step (automatically detected if present).
* --map : use maximum a-posteriori fitting (may significantly increase runtime).
* --cocluster : cluster SVs and SNVs together.
* --no_adjust : do not adjust read counts based on different classes of SV events.
* --subsample <integer> : (SNVs only); subsample N variants from total filtered input.
* --ss_seed (only valid if using subsampling) integer seed, can be set to ensure subsampling is replicable.
* --seeds <one integer per run, comma separated> : set a random seed per run in order to replicate pyMC runs.
* --XX and --XY : overwrite the config file genotype with XX or XY.

### Post-assigning SVs ###

The post-assign step obtains all the SVs or SNVs which were not used in the clustering (filtered out or not present due to sub-sampling) and assigns each variant to its most likely cluster, based on its read count and copy-number state. The step takes all the same input as the filter step (note: the -i flag is replaced with --svs). Optionally, the step also takes the --XX and XY config-override parameters specified in the cluster step. Additionally, the following parameters should also be defined if the <sample>_filtered_[snvs/svs].tsv file differs from the detault name and location (for instance, if sub-sampling was used):

* --filt_svs : filtered SV file that was used to run clustering. Defaults to <out>/<sample>_filtered_svs.tsv
* --filt_snvs : filtered SNV file that was used to run clustering. Defaults to <out>/<sample>_filtered_snvs.tsv

Additionally, the run directory may be specified:

* --run <run_dir> or -r <run_dir> : run for which to perform post-assignment. Defaults to "best" (selects best run if using MAP), if no best run is available, defaults to run0.

#### Post-assign output ####

Post-assign creates the same output structure as in the cluster step. One minor difference: in the assignment_probabilities file output, we use a transformed likelihood to estimate probability of a variant belonging to a particular cluster, rather than a probability based on the MCMC chain (which does not exist for a post-assigned variant).

#### Configuration file ####

If customisation is required for parameters, these can be modified in the svclone_config.ini, or a new config file can be specified and each step run with the --config or -cfg flag for all steps.

#### Post-processing (visualisation) ####

A post-processing script is included to visualise output of individual run plots and a summary plot for all runs of a sample. Run the script as follows:

    Rscript post_processing_fit_diagnostics.R <SVclone output dir> <sample name> <cnvs> [--map] [--snvs] [--coclus]

Optional inputs:

cnvs: for plotting circos plots (will draw copy-number tracks).

Optional flags:

* --map : use if SVclone was run with the map option (recommended if running multiple runs). This option will add fit metric plots.
* --snvs : use if SVclone output is for SNVs rather than SVs.
* --coclus: use if SNVs and SVs were coclustered.
