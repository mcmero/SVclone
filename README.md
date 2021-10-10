<p align="left">
<img src=https://github.com/mcmero/SVclone/blob/master/img/svclone_logo.png height=120/>
</p>

This package is used to cluster structural variants of similar cancer cell fraction (CCF). SVclone is divided into five components: annotate, count, filter, cluster and post-assign. The annotate step infers directionality of each breakpoint (if not supplied), recalibrates breakpoint position to the soft-clip boundary and subsequently classifies SVs using a rule-based approach. The count step counts the variant and non-variant reads from breakpoint locations. Both the annotate and count steps utilise BAM-level information. The filter step removes SVs based on a number of adjustable parameters and prepares the variants for clustering. SNVs can also be added at this step as well as CNV information, which is matched to SV and SNV loci. The post-assign step then allows SVs to be assigned to the derived model from SNV clustering. Optionally, SVs and SNVs can also be post-assigned to a joint SV + SNV model.

### How do I get set up? ###

First install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#). SVclone can then be installed via:

    conda install svclone -c bioconda -c conda-forge
    svclone --help

Alternatively, you may wish to install SVclone in its own conda virtual environment:

    conda create -n svclone -c bioconda -c conda-forge svclone
    conda activate svclone
    svclone --help

### Example data ###

Example data is provided to test your SVclone installation (data contains simulated clonal deletions). If you would like to run the tests, this can be done via:

    git clone https://github.com/mcmero/SVclone.git
    cd SVclone

    ./run_example.sh

You can check the following output plot, which will summarise the clustering result:

* tumour_p80_DEL/ccube_out/tumour

You can also test the simulated SNV data by running:

    ./run_example_wsnvs.sh

The simulated data contains a 100% CCF clone and a 30% subclone.

### Annotate step ###

Before running SVclone on real data, first download the configuration file via:

```
wget https://raw.githubusercontent.com/mcmero/SVclone/master/svclone_config.ini
```

Check the settings carefully and make sure that the config is approproate for your data set. Make sure that this file is in the directory from which you're running SVclone, or that you've specified the location using `-cfg <config_file>`. 

An indexed whole-genome sequencing BAM and a list of paired breakpoints from an SV caller of choice is required. This step is required for clustering of SVs, however, classifiation and directionality information from your choice of SV caller can be used rather than being inferred.

    svclone annotate -i <sv_input> -b <indexed_bamfile> -s <sample_name>

Input is expected in VCF format (directionality inferred from the ALT field is also supported). Each defined SV must have a matching mate, given in the MATEID value in the INFO section. Please see the [VCF spec](https://samtools.github.io/hts-specs/VCFv4.2.pdf) (section 3) for representing SVs using the VCF format. SVclone does not support unpaired break-ends, which means that the INFO field PARID must be specified (please see Section 5.4.4 in the VCF spec for an example).

Input may also be entered in Socrates or simple format (must be specified with `--sv_format simple` or `--sv_format socrates`). Simple format is as follows:

```
chr1	pos1	chr2	pos2
22	18240676	22	18232335
22	19940482	22	19937820
22	21383572	22	21382745
22	21395573	22	21395746
```

We recommend that directions from the SV caller of choice be used (`use_dir` must be set to `True` in the configuration file in this case). Optionally, if you already know the SV classifications, the name of the classification field can be specified in the config file (e.g. `sv_class_field: classification`).

An example of the 'full' SV simple format is as follows:

```
chr1	pos1	dir1	chr2	pos2	dir2	classification
22	18240676		-	22	18232335	-	INV
22	19940482		-	22	19937820    +	DEL
22	21383572		-	22	21382745	+	DUP
22	21395573		+	22	21395746	+	INV
```

Note that your classifications in your SV input will have to match those specified in the configuration file (these may be comma-separated):

```
[SVclasses]
# Naming conventions used to label SV types.
inversion_class: INV
deletion_class: DEL
dna_gain_class: DUP,INTDUP
dna_loss_class: DEL,INV,TRX
itrx_class: INTRX
```

Note that `dna_gain_class` will include any SV classification involving DNA duplication and `dna_loss_class` is any intra-chromosomal rearrangement *not involving* a gain (including balanced rearrangements). `itrx_class` refers to all inter-chromosomal translocations.

A blacklist (bed file) can also be supplied at this step to not process areas to remove SVs where any of its breakpoints fall into one of these areas.

#### Annotate Output ####

The above input example also corresponds with the output of this step (output to \<out\>/\<sample\>_svin.txt), with an added SV ID present in the column. Events that are considered part of the same event will have the same ID (which may be multiple breakpoints).

#### Required Parameters ####

* -i or --input : structural variants input file (see above).
* -b or --bam : bam file with corresponding index file.
* -s or --sample : Sample name. Will create processed output file as \<out\>/\<sample\>_svinfo.txt, parameters output as \<out\>/\<sample\>_params.txt.

#### Optional Parameters ####

* -o or --out \<directory\> : output directory to create files. Default: the sample name.
* -cgf or --config \<config.ini\> : SVclone configuration file with additional parameters (svclone_config.ini is the default).
* --sv_format \<vcf, simple, socrates\> : input format of SV calls, VCF by default, but may also be simple (see above) or from the SV caller Socrates.
* --blacklist \<file.bed\> : Takes a list of intervals in BED format. Skip processing of any break-pairs where either SV break-end overlaps an interval specified in the supplied bed file. Using something like the [DAC blacklist](https://www.encodeproject.org/annotations/ENCSR636HFF/) is recommended.

### Count step ###

Run SV processing submodule to obtain read counts around breakpoints on each sample BAM file like so:

    svclone count -i <svs> -b <indexed_bamfile> -s <sample_name>

The classification strings are not used by the program, except for DNA-gain events (such as duplications). The classification names for these types of SVs should be specified in the svclone_config.ini file (see configuration file section).

#### Count output ####

The count step will create a tab-separated \<out\>/\<sample\>_svinfo.txt file containing count information. For example:

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
* -s or --sample : Sample name. Will create processed output file as <out>/<sample>_svinfo.txt, parameters output as \<out\>/\<sample\>_params.txt.

#### Optional Parameters ####

* -o or --out \<directory\> : output directory to create files. Default: the sample name.
* -cgf or --config \<config.ini\>: SVclone configuration file with additional parameters (svclone_config.ini is the default).

### Filter step (Filter SVs and/or SNVs and attach CNV states) ###

To filter the data obtained from the SV counting program and/or filter SNV data, can be done like so:

    svclone filter -i <sv_info.txt> -s <sample_name>

Note that read length and insert sizes used by the filter step are provided as outputs from the count step (\<out\>/read_params.txt), based on the first 50,000 sampled reads in the bam file.

#### Filter output ####

The filter step outputs the file \<out\>/\<sample\>_filtered_svs.tsv and/or \<out\>/\<sample\>_filtered_snvs.tsv depending on input. For SVs, the output is akin to the _svinfo.txt file format with added fields:

* norm_mean: average of norm1 and norm2
* gtype1/2: copy-number states at loci 1 and 2: "major, minor, CNV fraction" for example, "1,1,1.0". May be subclonal if battenberg input is supplied e.g. "1,1,0.7|2,1,0.3".
* gtype1/2_adusted: the next closest copy-number state of the given locus in the same format as above.
* adjusted_norm1/2: normal read count at the given locus adjusted for DNA-gains.
* adjusted_support: total adjusted supporting read count.
* adjusted_depth1/2: adjusted_norm + adjusted_support at the given locus.
* raw_mean_vaf: support / (mean(norm1, norm2) + support)
* adjusted_vaf1/2: adjusted_support / adjusted_depth1/2

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

* -s or --sample \<name\> : sample name, currently only a single sample is supported. WARNING: if clustering using mutect SNVs, the sample name must match the sample name in the vcf file.
* -i or --input \<svinfo.txt\> : sv info file from SV count step.

#### Optional Parameters ####

* -o or --out \<out\> : output directory to create files. Default: the sample name.
* -cgf or --config \<config.ini\>: SVclone configuration file with additional parameters (svclone_config.ini is the default).
* --params \<params.txt\> : Parameters file from processing step containing read information. If not supplied, the default search path is \<out\>/\<sample\>_params.txt'
* -c or --cnvs \<cnv file\> : Battenberg subclones.txt or ASCAT caveman csv file containing segmented copy-numbers for patient sample. If not supplied, will assume that all regions are copy-number neutral (Major = 1, Minor = 1).
* -p \<file\> or --purity_ploidy \<file\>: Tumour purity and ploidy in tab-separated text file. Purity must be between 0 and 1. Ploidy can be a floating point number. Column names must be labelled 'sample', 'purity' and 'ploidy' (without quotes). Row 2 must contain the sample name and purity and ploidy values respectively.
* -g or --germline \<germline_svinfo.txt\> : Germline SVs; will filter out any tumour SVs which have >1 supporting read in the germline. Expects the same input format as the sv info file (you can run sv_process on the tumour SVs against the germline bam file to obtain this).
* --snvs \<snv_file\> : SNVs in VCF format to (optionally) compare the clustering with SVs.
* --snv_format \<sanger, mutect, mutect_callstats, consensus, multisnv\> (default = sanger) : Specify VCF input format (only if clustering SNVs).
* --blacklist \<file.bed\> : Takes a list of intervals in BED format. Skip processing of any break-pairs where either SV break-end overlaps an interval specified in the supplied bed file. Using something like the [DAC blacklist](https://www.encodeproject.org/annotations/ENCSR636HFF/) is recommended.

### SNV input formats ###

SVclone supports a number of SNV input formats specified by the `--snv_format` flag:

* sanger (default): VCF output from the [CaVEMan](http://cancerit.github.io/CaVEMan/) variant caller (allele depths contained in `FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ` fields).
* mutect: VCF output from [MuTect](https://software.broadinstitute.org/cancer/cga/mutect) ([example](https://github.com/mcmero/SVclone/blob/master/example_data/tumour_p80_DEL_snvs.vcf)).
* mutect_callstats: tab-delimited `call_stats.txt` output from MuTect.
* consensus: VCF consensus call format from the [PCAWG project](https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel).
* multisnv: multi-sample VCF output format from [multiSNV](https://bitbucket.org/joseph07/multisnv/wiki/Home).

### Copy-number input formats ###

Two copy-number input formats are supported: [Battenberg](https://github.com/cancerit/cgpBattenberg) and [ASCAT](https://github.com/cancerit/ascatNgs). The following shows a Battenberg subclones.txt example:

```
X	chr	startpos	endpos	BAF	pval	LogR	ntot	nMaj1_A	nMin1_A	frac1_A	nMaj2_A	nMin2_A	frac2_A
1	1	54491	7958088	0.509836471	1	0.029104607	2.510592339	1	1	1	NA	NA	NA
2	1	7958594	15094665	0.519197319	6.23E-36	-0.179467225	1.74352809	1	0	0.12550915	1	1	0.87449085
3	1	15095331	39549452	0.519007705	3.72E-09	-0.050483364	2.037000342	1	1	0.86585984	2	1	0.13414016
4	1	39556288	40141768	0.840924898	1	1.617870568	9.512437457	8	1	1	NA	NA	NA
```

Battenberg provides more possible segmentations, but SVclone will only use fractions for the '\_A' segments. Only the genomic coordinates, allelic copy-number and fraction fields are used by SVclone.

ASCAT's cavemam csv format is as follows:

```
1,1,54491,7958088,2,1,2,1
2,1,7958594,15094665,2,1,2,1
3,1,15095331,39549452,2,1,2,1
4,1,39556288,40141768,2,1,8,1
```

The fields are ID, chrom, start, end, normal total copy-number, normal minor allele copy-number, tumour total copy-number and tumour minor copy-number. This format gives copy-number calls at clonal resolution only.

### Purity/ploidy file ###

Purity/ploidy file must have the following layout (tab-separated):

```
sample	purity	ploidy
<sample>	<float between 0 - 1>	<positive float or integer>
```
### Clustering SVs ###

Once we have the filtered SV and/or SNV counts, we can run the clustering:

    svclone cluster -s <sample_name>

NOTE: cluster step must be run from the SVclone base directory.

#### Cluster output ####

SVclone creates output based on the PCAWG output specification. This includes (per run):

* number_of_clusters: number of clusters found.
* \<sample\>_multiplicity.txt: the total copy-number, the number of copies the variant occurs on, the different multiplicity options and the probability of each.
* \<sample\>_assignment_probability_table.txt: probability of each variant's assignment to each cluster, based on number of times the proportion that a variant occurs in a particular cluster over all MCMC iterations.
* \<sample\>_cluster_certainty.txt: each variant's most likely assignment to a particular cluster and corresponding average proportion (CCF x purity).
* \<sample\>_subclonal_structure: clusters found, the number of variants per cluster, the proportion and CCF.

Ccube will create an RData file under \<ccube_out\>/\<sample\>_[sv/snv]_results.RData

#### Required Parameters ####

* -s or --samples : sample names, currently only a single sample is supported.

#### Optional Parameters ####

* -o or --out : output directory (sample name by default)
* -cgf or --config \<config.ini\>: SVclone configuration file with additional parameters (svclone_config.ini is the default).
* --params \<params.txt\> : Parameters file from processing step containing read information. If not supplied, the default search path is \<out\>/\<sample\>_params.txt'
* --snvs \<filtered_snvs\> : SNVs output from filter step (automatically detected if present).
* --no_adjust : do not adjust read counts based on different classes of SV events.
* --subsample \<integer\> : (SNVs only); subsample N variants from total filtered input.
* --ss_seed (only valid if using subsampling) integer seed, can be set to ensure subsampling is replicable.
* --XX and --XY : overwrite the config file genotype with XX or XY (useful for bulk runs).

### Post-assigning SVs ###

Post-assignment involves reassigning SV cluster memberships to either an SNV model, or to a joint SV + SNV model. If using the joint-model approach, SNVs may also be reassigned using the derived joint model. Run this step as follows:

    svclone postassign -s <sample> --svs <sv_results.RData> --snvs <snv_results.RData> --out <output_dir> [--joint]

If the `--joint` argument is used, SVs and SNVs will be reassigned to a joint model. If this flag is omitted, SVs will be assigned to the SNV model (default behaviour).

#### Post-assign output ####

Post-assign creates the same output as the cluster step.

#### Configuration file ####

If customisation is required for parameters, these can be modified in the svclone_config.ini, or a new config file can be specified and each step run with the --config or -cfg flag for all steps.

### Data Availability (SVclone paper) ###

In silico sample mixtures were generated from patient data derived from patient 001 from the Hong et al. study (doi:10.1038/ncomms7605). The data are available in the EGA Sequence Read Archive under accession [EGAS00001000942](https://www.ebi.ac.uk/ega/studies/EGAS00001000942). Somatic and germline variant calls, mutational signatures, subclonal reconstructions, transcript abundance, splice calls and other core data generated by the ICGC/TCGA Pan-cancer Analysis of Whole Genomes Consortium is described here and available for download at [https://dcc.icgc.org/releases/PCAWG](https://dcc.icgc.org/releases/PCAWG). Additional information on accessing the data, including raw read files, can be found at https://docs.icgc.org/pcawg/data/. In accordance with the data access policies of the ICGC and TCGA projects, most molecular, clinical and specimen data are in an open tier which does not require access approval. To access controlled PCAWG data, please see the [ICGC DACO](https://daco.icgc.org/). Please see the following links for available PCAWG data related to the SVclone paper:

* [Consensus SVs](https://dcc.icgc.org/releases/PCAWG/consensus_sv) (open)
* [Consensus SNVs and INDELs](https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel) (open)
* [Consensus copy-numbers](https://dcc.icgc.org/releases/PCAWG/consensus_cnv) (open)
* [Subclonal reconstruction](https://dcc.icgc.org/releases/PCAWG/subclonal_reconstruction) (controlledl; includes SV CCFs)

All the other data supporting the findings of this study are available within the article and its supplementary information files and from the corresponding author upon reasonable request. A reporting summary for this article is available as a Supplementary Information file.

### Other code (SVclone paper) ###

Ccube clustering code can be found under [https://github.com/keyuan/ccube](https://github.com/keyuan/ccube). Code for generating all figures in the manuscript and the in silico mixture samples can be found under [https://github.com/mcmero/SVclone_Rmarkdown](https://github.com/mcmero/SVclone_Rmarkdown). Code for simulating SVs can be found under [https://github.com/mcmero/sv_simu_pipe](https://github.com/mcmero/sv_simu_pipe).
The core computational pipelines used by the PCAWG Consortium for alignment, quality control and variant calling are available to the public at [https://dockstore.org/search?search=pcawg](https://dockstore.org/search?search=pcawg) under the GNU General Public License v3.0, which allows for reuse and distribution.

#### NOTE: The master branch implements clustering using the [ccube](https://www.biorxiv.org/content/10.1101/484402v1.abstract) clustering model. If you would like to run the PyMC2 clustering model (described in the [biorXiv version](https://www.biorxiv.org/content/10.1101/172486v1.abstract)) please use the `pymc2` branch, although we highly recommend running the ccube version. ####
