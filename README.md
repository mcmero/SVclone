# README #

The script obtains all information necessary from a bam file and a list of structural variations in order to accurately count variant allele frequencies (VAF) and detect subclones from SV data. 

### How do I get set up? ###

Ensure you have the following dependencies installed:

* [PySam](http://pysam.readthedocs.org/en/latest/)
* [Pandas](http://pandas.pydata.org/)
* [Pandasql](https://pypi.python.org/pypi/pandasql)
* [sqlite3](https://docs.python.org/2/library/sqlite3.html)
* [ipdb](https://pypi.python.org/pypi/ipdb)

If you want to run unit tests:

* [nose](https://nose.readthedocs.org/en/latest/)

Install like so:

    python setup.py install

Run like so:

    python -m proc_svs.cmd <svs.txt> <indexed bamfile> <header.cfg> <output name> <average coverage>

The structural variation input file must be in the following tab-separated format:

```
bp1_chr	bp1_pos	bp1_dir	bp2_chr	bp2_pos	bp2_dir	classification
22	18240676	-	22	18232335	-	INV
22	19940482	-	22	19937820	-	INV
22	21383572	+	22	21382745	+	INV
22	21383573	-	22	21382746	-	INV 
```

The column names can be different, but must be specified in the header.cfg file, which looks like:

```
bp1_chr=bp1_chr
bp1_pos=bp1_pos
bp1_dir=bp1_dir
bp2_chr=bp2_chr
bp2_pos=bp2_pos
bp2_dir=bp2_dir
classification=classification
```

The left fields are used by the program (do not change these), the right fields correspond to the equivalent column in the SV text file.

### Who do I talk to? ###

For queries, contact cmerom@gmail.com
